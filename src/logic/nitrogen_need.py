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

# Directional kinetics for Arabidopsis GS isoforms. Vmax is relative
# (GLN1-2 normalised to 1.0). Compartment tracks the cellular location,
# matching the ``compartment`` field in the parser's metadata header.
# References:
#   Ishiyama et al. 2004, J Biol Chem 279:16598 -- GS1 isoform kinetics
#   Konishi et al. 2017, Plant J 91:913 -- GLN1-2 / 1-3 expression
#   Taira et al. 2004, Plant Cell 16:2048 -- GS2 Km (photorespiratory)
ISOFORM_KINETICS: dict[str, dict] = {
    "GLN1-2": {"km_nh4_uM": 260.0, "vmax_rel": 1.00, "compartment": "cytosol"},
    "GLN1-3": {"km_nh4_uM":  38.0, "vmax_rel": 0.35, "compartment": "cytosol"},
    "GLN2":   {"km_nh4_uM":  30.0, "vmax_rel": 3.20, "compartment": "plastid"},
}

# Two flux streams, one per compartment:
#
#   Cytosolic (GS1): primary ammonium assimilation from roots. Directly
#   responds to rhizosphere NH4+ -- this is what a fertiliser application
#   moves.
#
#   Plastidic (GS2): photorespiratory ammonia reassimilation. Driven by
#   light, temperature, and CO2 (which set the photorespiratory load),
#   not by soil N. Fertiliser application does not help here; if GS2 is
#   the bottleneck, more N in soil just increases losses.
#
# The recommendation target is the CYTOSOLIC mean rate -- the flux that
# responds to what the farmer can control.
TARGET_CYTOSOL_RATE = 0.80
DOSE_PER_DEFICIT_KG_HA = 40.0  # kg N/ha per unit of rate deficit


@dataclass
class SoilTelemetry:
    nh4_uM: float
    no3_uM: float
    moisture_pct: float
    ph: float
    temperature_c: float
    # Proxies for photorespiratory load (what GS2 sees). Photorespiration
    # scales with light * temperature / CO2; we pass rough scalars so the
    # logic can synthesise a plastid-facing NH4+ concentration without
    # modelling photosynthesis itself.
    ppfd_umol_m2_s: float = 800.0   # photosynthetic photon flux density
    leaf_temp_c: float = 22.0

    @classmethod
    def dummy(cls, regime: str = "limited") -> "SoilTelemetry":
        if regime == "replete":
            return cls(nh4_uM=800.0, no3_uM=3000.0,
                       moisture_pct=42.0, ph=6.2, temperature_c=21.0,
                       ppfd_umol_m2_s=900.0, leaf_temp_c=23.0)
        if regime == "limited":
            return cls(nh4_uM=30.0, no3_uM=120.0,
                       moisture_pct=28.0, ph=6.5, temperature_c=19.0,
                       ppfd_umol_m2_s=700.0, leaf_temp_c=21.0)
        if regime == "intermediate":
            return cls(nh4_uM=180.0, no3_uM=600.0,
                       moisture_pct=35.0, ph=6.3, temperature_c=20.0,
                       ppfd_umol_m2_s=800.0, leaf_temp_c=22.0)
        if regime == "heat_stress":
            # Hot/bright day: photorespiration spikes -> GS2 flux spikes.
            # Soil N is fine; the limit shifts inside the plastid.
            return cls(nh4_uM=300.0, no3_uM=1200.0,
                       moisture_pct=32.0, ph=6.4, temperature_c=30.0,
                       ppfd_umol_m2_s=1800.0, leaf_temp_c=34.0)
        raise ValueError(f"unknown regime: {regime}")

    def photorespiratory_nh4_uM(self) -> float:
        """Proxy NH4+ the plastid sees from photorespiratory glycine decarb.

        Empirical scaling, not mechanistic: flux goes up with light and
        with temperature (both push the RuBisCO oxygenation rate). The
        absolute calibration is placeholder; the shape is what matters
        for the recommender's regime detection.
        """
        light_factor = self.ppfd_umol_m2_s / 800.0
        # Rubisco oxygenation roughly doubles every 10 C above 20 C.
        temp_factor = 2.0 ** ((self.leaf_temp_c - 20.0) / 10.0)
        return 40.0 * light_factor * temp_factor


@dataclass
class PrecisionApplication:
    regime: str
    isoform_rates: dict[str, float]
    compartment_rates: dict[str, float]  # {"cytosol": float, "plastid": float}
    dominant_isoform: str
    limiting_compartment: str
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
    structures we've characterised and we pin the recommendation to
    their SOP Instance UIDs.

    Two-compartment model:
      * Cytosolic isoforms see the rhizosphere NH4+ signal and drive
        the recommended dose.
      * Plastidic isoforms see a photorespiratory NH4+ proxy from
        light/temperature; they do NOT drive the dose (fertiliser
        doesn't reach the plastid lumen) but they can flag regimes
        where soil N is fine yet leaf assimilation is limited by
        photorespiratory load.
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

    nh4_plastid = telemetry.photorespiratory_nh4_uM()

    rates: dict[str, float] = {}
    for iso, p in gated.items():
        substrate = nh4_plastid if p["compartment"] == "plastid" else telemetry.nh4_uM
        rates[iso] = _mm_rate(substrate, p["km_nh4_uM"], p["vmax_rel"])

    # Compartment-level means.
    cytosol_rates = [r for iso, r in rates.items()
                     if gated[iso]["compartment"] == "cytosol"]
    plastid_rates = [r for iso, r in rates.items()
                     if gated[iso]["compartment"] == "plastid"]
    compartment_rates = {
        "cytosol": round(sum(cytosol_rates) / len(cytosol_rates), 4)
                   if cytosol_rates else 0.0,
        "plastid": round(sum(plastid_rates) / len(plastid_rates), 4)
                   if plastid_rates else 0.0,
    }

    dominant = max(rates, key=rates.get)

    # GS1 (cytosol) and GS2 (plastid) act on DIFFERENT N sources --
    # primary root uptake vs photorespiratory recycling -- so they are
    # parallel pathways rather than a serial cascade. The "limiting
    # compartment" is therefore not a bottleneck in the steady-state
    # sense; it's a diagnostic for *which* source the plant is most
    # reliant on right now.
    def _utilisation(compartment: str) -> float:
        isos = [iso for iso in rates if gated[iso]["compartment"] == compartment]
        if not isos:
            return -1.0
        actual = sum(rates[iso] for iso in isos)
        ceiling = sum(gated[iso]["vmax_rel"] for iso in isos)
        return actual / ceiling if ceiling else -1.0

    utilisations = {"cytosol": _utilisation("cytosol"),
                    "plastid": _utilisation("plastid")}
    # Report the compartment running nearer its Vmax as the dominant
    # flux stream at this moment -- useful diagnostic, not a dose gate.
    limiting = max(utilisations, key=utilisations.get)

    # Dose: purely a function of cytosolic deficit below target. Soil
    # N application only moves the cytosolic pool, so that's the only
    # lever the recommender can pull.
    if cytosol_rates:
        cyto_util = utilisations["cytosol"]
        util_deficit = max(0.0, TARGET_CYTOSOL_RATE - cyto_util)
        dose = round(DOSE_PER_DEFICIT_KG_HA * util_deficit, 2)
    else:
        dose = 0.0

    if telemetry.nh4_uM < 50:
        regime = "N-limited"
        rationale = (
            "Low rhizosphere NH4+: high-affinity GLN1-3 carries cytosolic "
            "assimilation. Recommend a small split dose -- enough to lift "
            "GLN1-2 off its high-Km floor without overshooting."
        )
    elif telemetry.nh4_uM > 500:
        regime = "N-replete"
        rationale = (
            "Replete NH4+: GLN1-2 bulk assimilation dominates the cytosolic "
            "pool. Hold or taper application to avoid runoff."
        )
    else:
        regime = "N-intermediate"
        rationale = (
            "Intermediate NH4+: both cytosolic isoforms contribute. Apply "
            "a moderate maintenance dose tracking the cytosolic deficit."
        )

    # Plastid diagnostic: flag when leaf-level photorespiratory
    # recycling is running hot. Doesn't change the dose (parallel
    # pathway) but tells the operator the plant is under photo-
    # respiratory load and to watch for light/heat stress signals.
    if plastid_rates and utilisations["plastid"] > 0.80:
        rationale += (
            " PLASTID NOTE: photorespiratory GS2 is running at "
            f"{utilisations['plastid']:.0%} of Vmax -- high light/temp is "
            "driving heavy photorespiratory recycling in the chloroplasts. "
            "Leaf-level NUE is already near ceiling; more soil N aids "
            "cytosolic assimilation but will not further raise leaf throughput."
        )

    return PrecisionApplication(
        regime=regime,
        isoform_rates={k: round(v, 4) for k, v in rates.items()},
        compartment_rates=compartment_rates,
        dominant_isoform=dominant,
        limiting_compartment=limiting,
        soil_n_available_uM=telemetry.nh4_uM + telemetry.no3_uM,
        recommended_dose_kg_per_ha=dose,
        rationale=rationale,
        referenced_sop_instances=sop_refs,
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--regime",
        choices=["replete", "limited", "intermediate", "heat_stress"],
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
