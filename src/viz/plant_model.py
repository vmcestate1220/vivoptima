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
    "GLN1-2": (0.15, 0.75, 0.25),   # green
    "GLN1-3": (0.15, 0.35, 0.85),   # blue
}

# Canonical vascular nodes per isoform, derived from known expression
# patterns: GLN1-2 is enriched in mature leaves and phloem companion
# cells; GLN1-3 is enriched in young leaves / shoot apex under N stress.
# Coordinates are scene-space (arbitrary, stem runs up +y).
NODE_POSITIONS: dict[str, list[tuple[float, float, float]]] = {
    "GLN1-2": [(0.0, 1.0, 0.0), (0.0, 2.0, 0.0)],
    "GLN1-3": [(0.0, 3.0, 0.0), (0.0, 3.8, 0.0)],
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
        self.leaves = []
        for y, size in levels:
            for k, angle in enumerate(np.linspace(0, 2 * np.pi, 4, endpoint=False)):
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


def build_demo(
    pdb_dir: str | Path = "data/pdb/arabidopsis",
) -> PlantModel:
    """Parse every AlphaFold PDB in ``pdb_dir`` and attach it to the model."""
    pdb_dir = Path(pdb_dir)
    model = PlantModel()
    for pdb in sorted(pdb_dir.glob("AF-*.pdb")):
        parser = GSMonomerParser(pdb)
        _coords, ca_trace, header = parser.parse()
        nodes = NODE_POSITIONS.get(header.gene_symbol, [(0.0, 2.0, 0.0)])
        for node in nodes:
            model.attach_protein(header, ca_trace, node)
    return model


if __name__ == "__main__":
    build_demo()
    # VPython's local http server keeps the scene alive only while the
    # interpreter is running; this blocks so the browser tab stays up.
    input("Scene rendered. Press Enter to exit...")
