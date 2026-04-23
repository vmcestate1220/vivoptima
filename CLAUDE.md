# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository state

This is a **skeleton** repository. At the time of writing every source file (`src/main.py`, `docs/astm/standard_draft_v1.md`) is a 0-byte placeholder, `tests/` and `data/{structures,telemetry}/` are empty, and there is no `README`, `pyproject.toml`, `requirements.txt`, `Makefile`, or test runner configured. Before assuming anything, re-check file sizes — the scaffold is likely being filled in incrementally.

## Environment

A Python 3.13.2 virtual environment is checked into the repo at `./vivoptima/` (note: the venv directory is named the same as the project). Its own `.gitignore` excludes its contents from VCS.

Activate it before running any Python:

```bash
source vivoptima/bin/activate
# deactivate when done
```

Because there is no `requirements.txt` / `pyproject.toml`, the venv is currently the **only source of truth** for dependencies. If you add a package, also plan how it will be reproducibly installed (e.g. add a lockfile) rather than relying on the committed venv alone.

Key libraries already installed in the venv (domain signal — use these to orient, don't add parallel alternatives without reason):

- **Biopython** (`Bio`, `BioSQL`) — sequence / structural biology.
- **pydicom** — DICOM medical-imaging I/O.
- **VPython** + `jupyterlab_vpython` — in-notebook interactive 3D visualization (this is the viz stack for `src/viz/`, not matplotlib; only `matplotlib_inline` shim is present, not matplotlib itself).
- **numpy**, **pandas** — numerics / tabular data.
- **JupyterLab** (full stack incl. `jupyter_server_proxy`, `ipywidgets`) — primary interactive surface.

Launch the notebook UI with:

```bash
jupyter lab
```

## Layout intent

The empty directories encode the intended architecture; preserve these boundaries when adding code:

- `src/main.py` — entry point (currently empty).
- `src/metadata/` — metadata handling (given pydicom + ASTM context, likely DICOM tags and ASTM-standard descriptors).
- `src/viz/` — 3D visualization; target VPython given what's installed.
- `data/structures/` — 3D structural inputs (molecular or anatomical).
- `data/telemetry/` — time-series / measurement streams.
- `docs/astm/` — ASTM standard drafts authored in this repo (`standard_draft_v1.md`). Treat these as a first-class deliverable, not throwaway notes.
- `tests/` — test suite (no framework chosen yet; pytest is not currently installed — propose it before writing tests).

## Build / lint / test

None configured yet. If you are asked to run tests, build, or lint, **flag that these are not set up** rather than guessing a command — the user will likely want to decide the tooling.
