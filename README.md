 # Vivoptima

  ## Linked Contents

  - Overview (#overview)
  - Setup (#setup)
  - Quick Start (#quick-start)
  - Usage (#usage)
  - Dataset (#dataset)
  - References (#references)
  - Author (#author)

  ## Overview

  Vivoptima is a 3D "Nitrogen Digital Twin" platform designed to address the global fertilizer runoff crisis through
  "Clinical Agronomy." By mapping AlphaFold 3D structures of Glutamine Synthetase (GS) monomers, specifically the cytosolic
  GLN1.2 and GLN1.3 isoenzymes, onto a virtual model of Arabidopsis thaliana, Vivoptima predicts nitrogen requirements with
  molecular precision.

  The project leverages medical imaging logic (pydicom) and real-time 3D visualization (vpython) to create a Digital Twin
  Metadata Standard. This standard aims to bridge the gap between proteomic research and field-scale agricultural resource
  management, providing a framework for ASTM International committee consideration.

  ## Setup

  Ensure you have Python 3.11+ installed. It is recommended to use the provided virtual environment configuration.

  # Clone the repository
  git clone https://github.com/your-username/vivoptima.git
  cd vivoptima

  # Initialize environment
  python3 -m venv venv
  source venv/bin/activate

  # Install dependencies
  pip install pydicom vpython biopython numpy pandas

  ## Quick Start

  To initialize the digital twin and render the Arabidopsis model with default GS2 mapping:

  python src/main.py --model arabidopsis --isoenzyme gln1.2

  ## Usage

  Vivoptima functions by fusing three distinct data streams:

  1. Structural Mapping: Parsing .pdb files to locate active sites within the GS monomer.
  2. Metadata Tagging: Utilizing pydicom to wrap biological data in a clinical-grade header (e.g., mapping "Cultivar" to
     PatientName).
  3. Visual Simulation: Rendering the plant's vascular architecture in vpython to visualize nitrogen "sink" and "source"
     dynamics.

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