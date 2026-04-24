"""Parse AlphaFold GS monomer PDBs and emit a DICOM-analog metadata header.

Clinical analogy: a CT ROI carries coordinate data bound to a stable
SOP Instance UID that identifies one specific acquisition. Here each
protein conformation plays the same role -- the parsed structure is
the imaged volume, the active-site heavy atoms are the ROI contour,
and the UIDs let downstream modules (viz, logic) reference one exact
conformation without re-parsing.

Scope: handles all three Arabidopsis cytosolic/plastidic GS isoforms
currently in data/pdb/arabidopsis/:
  * Q8LCE1 (GLN1-2, AT1G66200, 356 aa, cytosol)
  * Q9LVI8 (GLN1-3, AT3G17820, 354 aa, cytosol)
  * Q43127 (GLN2,   AT5G35630, 430 aa precursor, plastid)
Expected residue counts and active-site residue numbers are pinned in
``ISOFORM_REGISTRY``; the parser guards against the wrong file being
fed for a registered accession. GS2's active-site residues are shifted
+58 from their GS1 positions (49-aa transit peptide + 9-aa N-terminal
extension in the mature body), verified by conserved-motif alignment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import uuid
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Chain import Chain

# DICOM PS3.5 B.2 blesses "2.25.<uuid128-as-int>" as a UID derivation
# scheme, letting us mint reproducible UIDs without registering an
# enterprise OID arc.
UID_ROOT = "2.25."

# Stable namespace for Vivoptima-minted UUIDv5s. Fixed so that every run
# of the parser on the same inputs produces the same UIDs.
VIVOPTIMA_NS = uuid.UUID("b0c6f3e2-5a11-4d3e-9d50-0cafebabe0a0")

# Our analog to a DICOM CT Image SOP Class: every parsed protein
# conformation claims this class. Minted once, reused forever.
SOP_CLASS_PROTEIN_CONFORMATION = UID_ROOT + str(
    uuid.uuid5(VIVOPTIMA_NS, "vivoptima.ProteinConformation").int
)

# Metal-binding pocket of plant cytosolic GS1 (Unno et al. 2006,
# J Biol Chem 281:29287, PDB 2D3A, Zea mays GS1). Residue numbers refer
# to the mature chain, 1-indexed. Verified conserved at the same
# positions in both Q8LCE1 and Q9LVI8.
GS1_ACTIVE_SITE_DEFAULT: tuple[int, ...] = (131, 192, 249, 297)

# GS2 (Q43127) is the plastidic/chloroplastic isoform. AlphaFold models
# the full-length precursor including the 49-aa plastid transit peptide,
# so the active-site residues sit +58 from their GS1 positions (the +58
# offset was verified by motif alignment across all four conserved
# residues -- the extra 9 aa beyond the 49-aa transit reflects a
# slightly longer N-terminus in the mature GS2 body).
GS2_ACTIVE_SITE_DEFAULT: tuple[int, ...] = tuple(p + 58 for p in GS1_ACTIVE_SITE_DEFAULT)

# Compartment labels feed into the Frame of Reference UID derivation so
# cytosolic and plastidic isoforms land in distinct reference frames --
# the DICOM-sanctioned way to declare "these coordinates live in
# different anatomical spaces."
COMPARTMENT_CYTOSOL = "cytosol"
COMPARTMENT_PLASTID = "plastid"

# Registry of isoforms we know how to parse. Extend, don't edit.
ISOFORM_REGISTRY: dict[str, dict] = {
    "Q8LCE1": {
        "gene_symbol": "GLN1-2",
        "locus_tag": "AT1G66200",
        "active_site": GS1_ACTIVE_SITE_DEFAULT,
        "compartment": COMPARTMENT_CYTOSOL,
        "expected_residues": 356,
    },
    "Q9LVI8": {
        "gene_symbol": "GLN1-3",
        "locus_tag": "AT3G17820",
        "active_site": GS1_ACTIVE_SITE_DEFAULT,
        "compartment": COMPARTMENT_CYTOSOL,
        "expected_residues": 354,
    },
    "Q43127": {
        "gene_symbol": "GLN2",
        "locus_tag": "AT5G35630",
        "active_site": GS2_ACTIVE_SITE_DEFAULT,
        "compartment": COMPARTMENT_PLASTID,
        "expected_residues": 430,
        "transit_peptide_end": 49,  # UniProt feature TRANSIT 1..49
    },
}


def _make_uid(*parts: str) -> str:
    payload = "::".join(parts)
    return UID_ROOT + str(uuid.uuid5(VIVOPTIMA_NS, payload).int)


@dataclass(frozen=True)
class ActiveSiteResidue:
    residue_number: int
    residue_name: str
    ca_coord: tuple[float, float, float]
    # Flattened heavy-atom xyz, matching DICOM Contour Data (3006,0050)
    # which is also a flat list of triplets.
    contour_data: list[float]

    @classmethod
    def from_biopython(cls, residue) -> "ActiveSiteResidue":
        ca = residue["CA"].coord
        contour: list[float] = []
        for atom in residue:
            if atom.element == "H":
                continue
            contour.extend(float(c) for c in atom.coord)
        return cls(
            residue_number=int(residue.id[1]),
            residue_name=residue.get_resname(),
            ca_coord=(float(ca[0]), float(ca[1]), float(ca[2])),
            contour_data=contour,
        )


@dataclass
class MetadataHeader:
    """JSON-serializable DICOM-analog header for one protein conformation."""

    # Patient IE -- organism is the "patient"
    patient_id: str
    patient_name: str
    # Study IE -- locus-level grouping
    study_instance_uid: str
    study_description: str
    # Series IE -- prediction run
    series_instance_uid: str
    series_description: str
    modality: str
    # Instance IE -- this exact conformation
    sop_class_uid: str
    sop_instance_uid: str
    frame_of_reference_uid: str
    # Equipment IE
    manufacturer: str
    software_version: str
    # ROI analog -- active-site contour
    structure_set_label: str
    roi_name: str
    roi_interpreted_type: str
    active_site_residues: list[ActiveSiteResidue]
    # Volume-level summary
    residue_count: int
    atom_count: int
    bounding_box: dict
    coordinate_checksum: str
    # Cross-references (non-DICOM, but handy downstream)
    source_file: str
    uniprot_accession: str
    gene_symbol: str
    locus_tag: str
    compartment: str

    def to_dict(self) -> dict:
        return asdict(self)

    def to_dicom_tagged(self) -> dict:
        """Re-key the header under canonical DICOM (gggg,eeee) tags.

        Non-DICOM fields land in a ``_vivoptima`` side-channel so the
        output still round-trips cleanly through JSON.
        """
        return {
            "(0010,0020)": self.patient_id,             # Patient ID
            "(0010,0010)": self.patient_name,           # Patient's Name
            "(0020,000D)": self.study_instance_uid,     # Study Instance UID
            "(0008,1030)": self.study_description,      # Study Description
            "(0020,000E)": self.series_instance_uid,    # Series Instance UID
            "(0008,103E)": self.series_description,     # Series Description
            "(0008,0060)": self.modality,               # Modality
            "(0008,0016)": self.sop_class_uid,          # SOP Class UID
            "(0008,0018)": self.sop_instance_uid,       # SOP Instance UID
            "(0020,0052)": self.frame_of_reference_uid, # Frame of Reference UID
            "(0008,0070)": self.manufacturer,           # Manufacturer
            "(0018,1020)": self.software_version,       # Software Versions
            "(3006,0002)": self.structure_set_label,    # Structure Set Label
            "(3006,0026)": self.roi_name,               # ROI Name
            "(3006,00A4)": self.roi_interpreted_type,   # RT ROI Interpreted Type
            "(3006,0050)": [
                {"residue_number": r.residue_number,
                 "residue_name": r.residue_name,
                 "contour_data": r.contour_data}
                for r in self.active_site_residues
            ],
            "_vivoptima": {
                "uniprot_accession": self.uniprot_accession,
                "gene_symbol": self.gene_symbol,
                "locus_tag": self.locus_tag,
                "compartment": self.compartment,
                "residue_count": self.residue_count,
                "atom_count": self.atom_count,
                "bounding_box": self.bounding_box,
                "coordinate_checksum": self.coordinate_checksum,
                "source_file": self.source_file,
            },
        }


class GSMonomerParser:
    """Parse one AlphaFold monomer PDB into coords + metadata header."""

    def __init__(
        self,
        pdb_path: str | Path,
        *,
        uniprot: str | None = None,
        chain_id: str = "A",
        expected_residue_count: int | None = None,
        active_site_residues: tuple[int, ...] | None = None,
    ) -> None:
        self.pdb_path = Path(pdb_path)
        self.chain_id = chain_id
        self.expected_residue_count = expected_residue_count

        self.uniprot = uniprot or self._infer_uniprot(self.pdb_path)
        reg = ISOFORM_REGISTRY.get(self.uniprot, {})
        self.gene_symbol = reg.get("gene_symbol", "UNKNOWN")
        self.locus_tag = reg.get("locus_tag", "UNKNOWN")
        self.compartment = reg.get("compartment", COMPARTMENT_CYTOSOL)
        self.active_site_residues = tuple(
            active_site_residues
            or reg.get("active_site")
            or GS1_ACTIVE_SITE_DEFAULT
        )
        # If the caller didn't pin a residue count and the registry has
        # one, use the registry value so registered isoforms are guarded
        # automatically.
        if self.expected_residue_count is None and "expected_residues" in reg:
            self.expected_residue_count = reg["expected_residues"]

    @staticmethod
    def _infer_uniprot(path: Path) -> str:
        # AlphaFold filenames: AF-<accession>-F1-model_v*.pdb
        stem = path.stem
        if stem.startswith("AF-"):
            return stem.split("-")[1]
        raise ValueError(
            f"Cannot infer UniProt accession from {path.name}; "
            f"pass uniprot=... explicitly"
        )

    def parse(self) -> tuple[np.ndarray, np.ndarray, MetadataHeader]:
        """Return ``(heavy_atom_coords, ca_trace, header)``.

        ``ca_trace`` is shape ``(N_residues, 3)`` and is extracted here
        so downstream viz doesn't need to re-walk the PDB or guess a
        stride through mixed heavy atoms.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(self.uniprot, str(self.pdb_path))
        model = next(structure.get_models())
        chain: Chain = model[self.chain_id]

        residues = [r for r in chain if r.id[0] == " "]
        coords = np.array(
            [a.coord for r in residues for a in r if a.element != "H"],
            dtype=np.float64,
        )
        ca_trace = np.array(
            [r["CA"].coord for r in residues if "CA" in r],
            dtype=np.float64,
        )

        if self.expected_residue_count is not None \
                and len(residues) != self.expected_residue_count:
            raise ValueError(
                f"{self.pdb_path.name}: expected "
                f"{self.expected_residue_count} residues, found {len(residues)}"
            )

        by_resnum = {r.id[1]: r for r in residues}
        missing = [n for n in self.active_site_residues if n not in by_resnum]
        if missing:
            raise KeyError(
                f"{self.pdb_path.name}: active-site residues {missing} absent"
            )
        active = [
            ActiveSiteResidue.from_biopython(by_resnum[n])
            for n in self.active_site_residues
        ]

        checksum = hashlib.sha256(coords.tobytes()).hexdigest()
        bbox = {
            "min": coords.min(axis=0).tolist(),
            "max": coords.max(axis=0).tolist(),
        }

        header = MetadataHeader(
            patient_id="taxid:3702",
            patient_name="Arabidopsis thaliana",
            study_instance_uid=_make_uid(
                "Study", "ArabidopsisThaliana", "NitrogenAssimilation"
            ),
            study_description="Vivoptima nitrogen assimilation digital twin",
            series_instance_uid=_make_uid(
                "Series", "AlphaFold-v6", self.locus_tag
            ),
            series_description=(
                f"AlphaFold monomer prediction: {self.gene_symbol}"
            ),
            modality="AF",  # AlphaFold -- not a real DICOM modality code
            sop_class_uid=SOP_CLASS_PROTEIN_CONFORMATION,
            sop_instance_uid=_make_uid(
                "SOPInstance", self.uniprot, checksum
            ),
            frame_of_reference_uid=_make_uid(
                "FrameOfReference", self.compartment, self.uniprot, checksum[:16]
            ),
            manufacturer="DeepMind / EBI AlphaFold DB",
            software_version="AlphaFold Monomer v2.0 (model v6)",
            structure_set_label=f"{self.gene_symbol} catalytic pocket",
            roi_name=f"{self.gene_symbol} metal-binding pocket",
            roi_interpreted_type="CATALYTIC_SITE",
            active_site_residues=active,
            residue_count=len(residues),
            atom_count=int(coords.shape[0]),
            bounding_box=bbox,
            coordinate_checksum=checksum,
            source_file=str(self.pdb_path.resolve()),
            uniprot_accession=self.uniprot,
            gene_symbol=self.gene_symbol,
            locus_tag=self.locus_tag,
            compartment=self.compartment,
        )
        return coords, ca_trace, header


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb", nargs="+", type=Path, help="PDB file(s) to parse")
    ap.add_argument(
        "--out-dir", type=Path,
        default=Path("data/metadata"),
        help="Directory to write JSON headers into",
    )
    ap.add_argument(
        "--expected-residues", type=int, default=None,
        help="Assert residue count (e.g. 430 for GS2; omit for GS1 files)",
    )
    ap.add_argument(
        "--dicom-tagged", action="store_true",
        help="Emit keys as DICOM (gggg,eeee) notation",
    )
    args = ap.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for path in args.pdb:
        parser = GSMonomerParser(
            path, expected_residue_count=args.expected_residues
        )
        coords, _ca, header = parser.parse()
        payload = (
            header.to_dicom_tagged() if args.dicom_tagged else header.to_dict()
        )
        out = args.out_dir / f"{path.stem}.meta.json"
        out.write_text(json.dumps(payload, indent=2, default=str))
        print(
            f"{path.name}: {len(coords)} heavy atoms, "
            f"{header.residue_count} residues -> {out}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
