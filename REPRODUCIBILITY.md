# Reproducibility checklist

This repository intentionally does not track the full generated CSV datasets. They are large derived artifacts, and they can be regenerated from the deterministic scripts in this repository.

The commands below assume they are run from the repository root.

## Environment

Use Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
