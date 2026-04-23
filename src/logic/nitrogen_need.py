"""Precision-application recommender.

Compares soil telemetry against the parsed GS1 isoenzyme metadata
(from ``src.metadata.pdb_parser``) and emits a ``PrecisionApplication``
value containing a recommended N dose and per-isoform rate contribution.

Scope today: the kinetic parameters are literature-derived scalars; the
structural metadata is used to *gate* which isoforms are in scope (we
only make recommendations against isoforms whose conformations we've
parsed) and to pin each recommendation to specific SOP Instance UIDs
for reproducibility. A future pass can derive kinetics directly from
the active-site geometry in the metadata header.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

# Directional kinetics for Arabidopsis cytosolic GS1 isoforms. Vmax is
# relative (GLN1-2 normalised to 1.0). References:
#   Ishiyama et al. 2004, J Biol Chem 279:16598 -- GS1 isoform kinetics
#   Konishi et al. 2017, Plant J 91:913 -- GLN1-2/1-3 expression patterns
ISOFORM_KINETICS: dict[str, dict[str, float]] = {
    "GLN1-2": {"km_nh4_uM": 260.0, "vmax_rel": 1.00},
    "GLN1-3": {"km_nh4_uM":  38.0, "vmax_rel": 0.35},
}

# Agronomic threshold: aim for the *mean* isoform velocity to sit at or
# above 80% of its Vmax. Below, we recommend applying N; above, taper
# to minimise leaching runoff.
TARGET_MEAN_RATE = 0.80
DOSE_PER_DEFICIT_KG_HA = 40.0  # kg N/ha per unit of rate deficit


@dataclass
class SoilTelemetry:
    nh4_uM: float
    no3_uM: float
    moisture_pct: float
    ph: float
    temperature_c: float

    @classmethod
    def dummy(cls, regime: str = "limited") -> "SoilTelemetry":
        if regime == "replete":
            return cls(nh4_uM=800.0, no3_uM=3000.0,
                       moisture_pct=42.0, ph=6.2, temperature_c=21.0)
        if regime == "limited":
            return cls(nh4_uM=30.0, no3_uM=120.0,
                       moisture_pct=28.0, ph=6.5, temperature_c=19.0)
        if regime == "intermediate":
            return cls(nh4_uM=180.0, no3_uM=600.0,
                       moisture_pct=35.0, ph=6.3, temperature_c=20.0)
        raise ValueError(f"unknown regime: {regime}")


@dataclass
class PrecisionApplication:
    regime: str
    isoform_rates: dict[str, float]
    dominant_isoform: str
    soil_n_available_uM: float
    recommended_dose_kg_per_ha: float
    rationale: str
    referenced_sop_instances: dict[str, str]


def _mm_rate(nh4_uM: float, km: float, vmax: float) -> float:
    return vmax * nh4_uM / (km + nh4_uM)


def _load_headers(metadata_dir: Path) -> dict[str, dict]:
    """Load every ``*.meta.json`` written by the parser, keyed by gene."""
    headers: dict[str, dict] = {}
    if not metadata_dir.exists():
        return headers
    for p in sorted(metadata_dir.glob("*.meta.json")):
        data = json.loads(p.read_text())
        # Parser emits either a flat header or a DICOM-tagged dict; tolerate both.
        gene = (
            data.get("gene_symbol")
            or data.get("_vivoptima", {}).get("gene_symbol")
        )
        if gene:
            headers[gene] = data
    return headers


def _sop_instance(header: dict) -> str:
    return (
        header.get("sop_instance_uid")
        or header.get("(0008,0018)")
        or "UNKNOWN"
    )


def recommend(
    telemetry: SoilTelemetry,
    isoform_headers: dict[str, dict] | None = None,
) -> PrecisionApplication:
    """Return a precision-application recommendation.

    ``isoform_headers`` maps gene symbol -> metadata header dict from
    the parser. When supplied, we only evaluate isoforms whose
    structures we've characterised (keeps the loop grounded in the
    actual conformations in ``data/pdb/arabidopsis/``) and we pin the
    recommendation to their SOP Instance UIDs.
    """
    if isoform_headers:
        gated = {g: k for g, k in ISOFORM_KINETICS.items() if g in isoform_headers}
        sop_refs = {g: _sop_instance(h) for g, h in isoform_headers.items()
                    if g in gated}
    else:
        gated = ISOFORM_KINETICS
        sop_refs = {}

    if not gated:
        raise ValueError(
            "No isoforms to evaluate; metadata headers and ISOFORM_KINETICS "
            "have no overlap"
        )

    rates = {
        iso: _mm_rate(telemetry.nh4_uM, p["km_nh4_uM"], p["vmax_rel"])
        for iso, p in gated.items()
    }
    dominant = max(rates, key=rates.get)
    mean_rate = sum(rates.values()) / len(rates)
    deficit = max(0.0, TARGET_MEAN_RATE - mean_rate)
    dose = round(DOSE_PER_DEFICIT_KG_HA * deficit, 2)

    if telemetry.nh4_uM < 50:
        regime = "N-limited"
        rationale = (
            "Low rhizosphere NH4+: high-affinity GLN1-3 carries assimilation. "
            "Recommend a small split dose -- enough to lift GLN1-2 off its "
            "high-Km floor without overshooting into leaching-risk territory."
        )
    elif telemetry.nh4_uM > 500:
        regime = "N-replete"
        rationale = (
            "Replete NH4+: GLN1-2 bulk assimilation dominates. Hold or taper "
            "application to avoid runoff and reactive-N loss."
        )
    else:
        regime = "N-intermediate"
        rationale = (
            "Intermediate NH4+: both isoforms contribute. Apply a moderate "
            "maintenance dose tracking the computed deficit."
        )

    return PrecisionApplication(
        regime=regime,
        isoform_rates={k: round(v, 4) for k, v in rates.items()},
        dominant_isoform=dominant,
        soil_n_available_uM=telemetry.nh4_uM + telemetry.no3_uM,
        recommended_dose_kg_per_ha=dose,
        rationale=rationale,
        referenced_sop_instances=sop_refs,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--regime",
        choices=["replete", "limited", "intermediate"],
        default="limited",
    )
    ap.add_argument(
        "--metadata-dir", type=Path, default=Path("data/metadata"),
        help="Where parser wrote *.meta.json headers",
    )
    args = ap.parse_args(argv)

    headers = _load_headers(args.metadata_dir)
    telemetry = SoilTelemetry.dummy(args.regime)
    rec = recommend(telemetry, headers or None)

    print(json.dumps(
        {
            "telemetry": asdict(telemetry),
            "recommendation": asdict(rec),
            "isoforms_in_scope": sorted(headers),
        },
        indent=2,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
