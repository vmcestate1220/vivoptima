"""Virtual Arabidopsis plant with GS1 isoform overlays in VPython.

Builds a simplified stem-and-leaves plant from geometric primitives and
attaches parsed GS1 monomers to canonical vascular nodes. Each isoform
is color-coded: GLN1-2 = green (bulk-flow, expressed in mature source
tissues), GLN1-3 = blue (high-affinity, induced under N limitation).

Runs in JupyterLab (inline via jupyterlab_vpython) or as a standalone
script (opens a browser canvas served by vpython's local http server).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from src.metadata.pdb_parser import GSMonomerParser, MetadataHeader

# Importing vpython is heavy and currently blocked by a setuptools 82
# / pkg_resources incompatibility; do it lazily so the module can still
# be imported (e.g. for type checks) in environments where vpython is
# not launchable.
def _vp():
    import vpython  # noqa: WPS433 (runtime import is intentional)
    return vpython


ISOFORM_COLORS_RGB: dict[str, tuple[float, float, float]] = {
    "GLN1-2": (0.15, 0.75, 0.25),   # green  -- cytosolic, mature leaves
    "GLN1-3": (0.15, 0.35, 0.85),   # blue   -- cytosolic, young tissue / stress
    "GLN2":   (0.95, 0.55, 0.10),   # amber  -- plastidic, photorespiration
}

# Canonical sites per isoform.
#   Cytosolic GS1 isoforms live along the vascular system (stem + phloem
#   companion cells), so they map to stem-axis nodes. GLN1-2 is enriched
#   in mature tissue (lower stem), GLN1-3 in young tissue / shoot apex.
#
#   Plastidic GS2 lives in chloroplasts inside mesophyll cells, so it
#   is rendered at multiple scattered points across each leaf surface
#   rather than at a single stem node. These are filled in at render
#   time from the leaf geometry (see PlantModel._leaf_chloroplast_sites).
VASCULAR_NODES: dict[str, list[tuple[float, float, float]]] = {
    "GLN1-2": [(0.0, 1.0, 0.0), (0.0, 2.0, 0.0)],
    "GLN1-3": [(0.0, 3.0, 0.0), (0.0, 3.8, 0.0)],
}

# How many chloroplast glyphs to scatter per leaf for GS2. One is
# enough to mark the isoform's presence in every leaf; more reads as
# visual noise because each glyph is a full 430-residue CA curve.
CHLOROPLASTS_PER_LEAF = 1

# Per-isoform render spec overrides. Plastidic glyphs sit *inside* the
# leaf ellipsoids (which are only ~0.04 m tall), so they need a much
# smaller scale and thinner backbone than vascular stem-node glyphs.
ISOFORM_SPEC_OVERRIDES: dict[str, dict] = {
    "GLN2": {
        "protein_scale": 0.006,
        "backbone_radius": 0.008,
        "active_site_radius": 0.025,
    },
}


@dataclass
class RenderSpec:
    protein_scale: float = 0.02
    protein_opacity: float = 0.55
    backbone_radius: float = 0.025
    active_site_radius: float = 0.09


class PlantModel:
    def __init__(self, title: str = "Vivoptima virtual Arabidopsis"):
        vp = _vp()
        self.vp = vp
        self.scene = vp.canvas(
            title=title, width=960, height=720,
            background=vp.vector(1, 1, 1),
        )
        self.scene.forward = vp.vector(-0.3, -0.2, -1)
        self.scene.center = vp.vector(0, 2.0, 0)
        self._build_stem()
        self._build_leaves()
        self._nodes_drawn: list[tuple[str, vp.vector]] = []

    def _build_stem(self) -> None:
        vp = self.vp
        self.stem = vp.cylinder(
            pos=vp.vector(0, 0, 0),
            axis=vp.vector(0, 4.2, 0),
            radius=0.08,
            color=vp.vector(0.25, 0.5, 0.2),
        )

    def _build_leaves(self) -> None:
        vp = self.vp
        # (y, radial_size) pairs, roughly phyllotactic.
        levels = [(0.8, 0.55), (1.8, 0.75), (2.8, 0.65), (3.6, 0.45)]
        self.leaves: list = []
        self._leaf_frames: list[tuple[tuple[float, float, float],
                                      tuple[float, float, float],
                                      float]] = []
        for y, size in levels:
            for angle in np.linspace(0, 2 * np.pi, 4, endpoint=False):
                direction = vp.vector(np.cos(angle), 0.0, np.sin(angle))
                mid = vp.vector(0, y, 0) + direction * (size * 0.5)
                self.leaves.append(
                    vp.ellipsoid(
                        pos=mid,
                        length=size, height=0.04, width=size * 0.4,
                        axis=direction,
                        color=vp.vector(0.35, 0.72, 0.32),
                    )
                )
                self._leaf_frames.append(
                    ((mid.x, mid.y, mid.z),
                     (direction.x, direction.y, direction.z),
                     size)
                )

    def _leaf_chloroplast_sites(self) -> list[tuple[float, float, float]]:
        """Scatter CHLOROPLASTS_PER_LEAF points across each leaf surface."""
        sites: list[tuple[float, float, float]] = []
        rng = np.random.default_rng(seed=42)  # deterministic placement
        for (mid, axis, size) in self._leaf_frames:
            mid_v = np.array(mid)
            axis_v = np.array(axis)
            # Build an in-plane perpendicular for width scatter
            perp = np.cross(axis_v, np.array([0.0, 1.0, 0.0]))
            if np.linalg.norm(perp) < 1e-6:
                perp = np.array([1.0, 0.0, 0.0])
            perp = perp / np.linalg.norm(perp)
            for _ in range(CHLOROPLASTS_PER_LEAF):
                u = rng.uniform(-0.4, 0.4) * size        # along leaf axis
                v = rng.uniform(-0.15, 0.15) * size      # across leaf
                p = mid_v + axis_v * u + perp * v
                sites.append((float(p[0]), float(p[1]), float(p[2])))
        return sites

    def attach_protein(
        self,
        header: MetadataHeader,
        ca_trace: np.ndarray,
        node_position: tuple[float, float, float],
        *,
        spec: RenderSpec = RenderSpec(),
    ) -> None:
        """Render one parsed monomer at ``node_position``.

        ``ca_trace`` is the ``(N_residues, 3)`` CA array from the
        parser. We stroke it as a VPython curve coloured by isoform
        and highlight active-site residues with larger emissive
        spheres in the same colour.
        """
        vp = self.vp
        rgb = ISOFORM_COLORS_RGB.get(header.gene_symbol, (0.5, 0.5, 0.5))
        col = vp.vector(*rgb)
        origin = vp.vector(*node_position)

        centroid = ca_trace.mean(axis=0)
        centered = (ca_trace - centroid) * spec.protein_scale

        trace = vp.curve(color=col, radius=spec.backbone_radius)
        for xyz in centered:
            trace.append(origin + vp.vector(*xyz))

        for res in header.active_site_residues:
            rel = (np.array(res.ca_coord) - centroid) * spec.protein_scale
            vp.sphere(
                pos=origin + vp.vector(*rel),
                radius=spec.active_site_radius,
                color=col,
                opacity=spec.protein_opacity + 0.3,
                emissive=True,
            )

        vp.label(
            pos=origin + vp.vector(0, 0.35, 0),
            text=f"{header.gene_symbol}\n{header.locus_tag}",
            height=12, color=col, box=False, opacity=0,
        )
        self._nodes_drawn.append((header.gene_symbol, origin))


def _target_nodes(model: PlantModel, header: MetadataHeader) \
        -> list[tuple[float, float, float]]:
    """Pick rendering positions for one isoform based on compartment.

    Cytosolic isoforms get the vascular stem nodes mapped by gene
    symbol. Plastidic isoforms get scattered across leaf surfaces
    (one glyph per chloroplast site).
    """
    if header.compartment == "plastid":
        return model._leaf_chloroplast_sites()
    return VASCULAR_NODES.get(header.gene_symbol, [(0.0, 2.0, 0.0)])


def build_demo(
    pdb_dir: str | Path = "data/pdb/arabidopsis",
) -> PlantModel:
    """Parse every AlphaFold PDB in ``pdb_dir`` and attach it to the model."""
    pdb_dir = Path(pdb_dir)
    model = PlantModel()
    for pdb in sorted(pdb_dir.glob("AF-*.pdb")):
        parser = GSMonomerParser(pdb)
        _coords, ca_trace, header = parser.parse()
        spec = RenderSpec(**ISOFORM_SPEC_OVERRIDES.get(header.gene_symbol, {}))
        for node in _target_nodes(model, header):
            model.attach_protein(header, ca_trace, node, spec=spec)
    return model


if __name__ == "__main__":
    build_demo()
    # VPython's local http server keeps the scene alive only while the
    # interpreter is running; this blocks so the browser tab stays up.
    input("Scene rendered. Press Enter to exit...")
