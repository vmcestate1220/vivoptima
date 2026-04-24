# Vivoptima

[![Python 3.13](https://img.shields.io/badge/python-3.13-blue.svg)](https://www.python.org/)
[![Biopython](https://img.shields.io/badge/Biopython-1.87-green.svg)](https://biopython.org/)
[![VPython](https://img.shields.io/badge/VPython-7.6-orange.svg)](https://vpython.org/)
[![AlphaFold DB](https://img.shields.io/badge/AlphaFold-DB-5555ff.svg)](https://alphafold.ebi.ac.uk/)
[![DICOM](https://img.shields.io/badge/DICOM-analog-9cf.svg)](https://www.dicomstandard.org/)
[![Status](https://img.shields.io/badge/status-early--scaffold-lightgrey.svg)]()

## Contents

- [Overview](#overview)
- [Setup](#setup)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Dataset](#dataset)
- [References](#references)
- [Author](#author)

## Overview

Vivoptima is a 3D "Nitrogen Digital Twin" platform designed to address the global fertilizer runoff crisis through "Clinical Agronomy." By mapping AlphaFold 3D structures of Glutamine Synthetase (GS) monomers, specifically the cytosolic GLN1.2 and GLN1.3 isoenzymes, onto a virtual model of Arabidopsis thaliana, Vivoptima predicts nitrogen requirements with molecular precision.

The project leverages medical imaging logic (pydicom) and real-time 3D visualization (vpython) to create a Digital Twin Metadata Standard. This standard aims to bridge the gap between proteomic research and field-scale agricultural resource management, providing a framework for ASTM International committee consideration.

## Setup

This project targets **Python 3.13**. A virtual environment ships at `./vivoptima/` (its contents are excluded from VCS by the venv's own `.gitignore`); most contributors will prefer a fresh venv built from `requirements.txt`.

```bash
git clone https://github.com/vmcestate1220/vivoptima.git
cd vivoptima

python3.13 -m venv .venv
source .venv/bin/activate

pip install -r requirements.txt
```

The `setuptools<81` pin in `requirements.txt` is required because VPython 7.6.5 imports `pkg_resources`, which setuptools removed in 81+. Drop the pin once VPython releases a version without that import.

## Quick Start

The three modules run as `python -m` entry points:

```bash
# 1. Parse the AlphaFold PDBs into DICOM-analog JSON metadata headers.
python -m src.metadata.pdb_parser data/pdb/arabidopsis/AF-*.pdb

# 2. Compute a precision-application recommendation from dummy soil telemetry.
python -m src.logic.nitrogen_need --regime limited

# 3. Open the 3D virtual plant in a browser (press Enter in the terminal to exit).
python -m src.viz.plant_model
```

Top-level orchestration is not wired up yet: `src/main.py` is still a placeholder.

## Usage

Vivoptima functions by fusing three distinct data streams:

1. Structural Mapping: Parsing .pdb files to locate active sites within the GS monomer.
2. Metadata Tagging: Utilizing pydicom to wrap biological data in a clinical-grade header (e.g., mapping "Cultivar" to PatientName).
3. Visual Simulation: Rendering the plant's vascular architecture in vpython to visualize nitrogen "sink" and "source" dynamics.

## Dataset

- Protein Structures: AlphaFold models for Arabidopsis thaliana (Locus AT1G66200 and AT3G17820).
- Imaging Data: Synthetic multispectral hypercubes (generated in data/synthetic) simulating NIR and Red-Edge reflectance.
- Standards: Draft metadata dictionaries located in docs/astm/.

## References

- Jumper, J. et al. (2021). "Highly accurate protein structure prediction with AlphaFold." Nature.
- DICOM Standard: WG-21 (Computed Tomography) for spatial metadata logic.
- ASTM Committee E62 on Industrial Biotechnology.
- Arabidopsis GS1 Isoenzyme Research: Studies on the distinct roles of GLN1.2 vs. GLN1.3 in nitrogen loading.

## Author

Christopher D. Cocchiaraley
