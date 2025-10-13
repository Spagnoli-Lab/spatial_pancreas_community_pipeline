# Spatial Pipeline Overview

This repository contains a Nextflow-driven workflow that takes raw spatial transcriptomics mask outputs and masked read coordinates, filters the segmentation, overlays QC plots, runs pciSeq to assign cell types, projects single-cell gene signatures with Tangram, and produces post-Tangram quality reports. The end-to-end run generates:

- Filtered segmentation masks and overlays for visual QC.
- pciSeq classification artifacts (pickled results and AnnData ready for downstream tools).
- Tangram-derived spatial AnnData (predicted labels, projected gene expression, mapping object).
- Tangram QC plots (training-score trends, annotation scatter plot) and a cell-type count report.

## Environment Setup

1. **Install Conda or Mamba** – the pipeline relies on per-process Conda environments. Mamba speeds up environment creation but is optional.
2. **Install Nextflow (>= 25.04.x)** – download the latest release or use `curl -s https://get.nextflow.io | bash` and move the binary into your `PATH`.
3. **Clone the repository**:
   ```bash
   git clone https://github.com/<your-org>/spatial_pipeline.git
   cd spatial_pipeline
   ```
4. **(Optional) Prepare input data** – update `samplesheet.csv` with paths to your DAPI masks, cellpose outputs, and masked read CSVs. Ensure the single-cell reference AnnData (`test_data/single_cell/Cartana_simplified_fixname.h5ad` by default) matches your experiment.

All other dependencies are installed automatically during the Nextflow run via the Conda environment specifications (`envs/*.yml`).

## Running the Pipeline

### Quick Test Run

Execute the pipeline using the bundled test data:
```bash
nextflow run main.nf -profile test
```

Key outputs are written under `outputs/nf/`:
- `cell_filter/` – filtered masks (`*_processed_mask.npz`, `_cell_masks.tif`).
- `overlay_qc/` – overlay PNGs and summaries.
- `pciseq/` – pciSeq results and Tangram inputs.
- `tangram/` – Tangram spatial AnnData, projected genes, mapping.
- `tangram/post/` – postprocessing plots and cell-type counts.

## Directory Structure

```
envs/                 # Conda environment specs per process
scripts/              # Python utilities (filtering, overlays, pciSeq, Tangram, plotting)
test_data/            # Demo inputs
samplesheet.csv       # Sample metadata for the test run
main.nf               # Nextflow workflow definition
nextflow.config       # Profiles and global settings
```

## Troubleshooting

- If a Conda environment fails to create, ensure Conda/Mamba is available in `PATH`.
- For plotting errors (e.g. missing GUI backend), the scripts save figures directly to disk—no display is required.
- Rerun with `-resume` after fixing issues to skip successful steps.

For details on each script, review the docstrings inside `scripts/*.py`.
