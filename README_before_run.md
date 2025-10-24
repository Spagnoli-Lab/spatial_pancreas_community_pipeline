# Before you run

This section describes how to prepare your input data and sample sheet so the pipeline can run without errors.

## 1) Expected directory structure (root = project folder)
- data/
  - dapi/                       -> DAPI images (DAPI_<sample>.tif)
  - masks/                      -> masks (DAPI_<sample>_mask.tif)
  - cellpose/                   -> segmentation arrays (DAPI_<sample>_seg.npy)
  - unmasked_reads/             -> per-sample reads CSV (<sample>.csv)
  - masked_reads/               -> per-sample masked reads CSV (<sample>_masked.csv)

Example (relative paths)
- data/dapi/DAPI_E12.5_9.tif
- data/masks/DAPI_E12.5_9_mask.tif
- data/cellpose/DAPI_E12.5_9_seg.npy
- data/unmasked_reads/E12.5_9.csv
- data/masked_reads/E12.5_9_masked.csv

## 2) File naming conventions
- sample id: use a simple token (no spaces). Example: E12.5_9
- DAPI image: DAPI_<sample>.tif
- Mask image: DAPI_<sample>_mask.tif
- Segmentation file: DAPI_<sample>_seg.npy
- Unmasked reads csv: <sample>.csv
- Masked reads csv: <sample>_masked.csv

Notes:
- Paths in the sample sheet can be relative to the project root or absolute; be consistent.
- Avoid spaces and special characters in filenames (use underscores).
- Ensure one segmentation file per sample and matching sample ids across files.

## 3) Preparing the sample sheet
- CSV with header and one row per sample.
- Required columns (order and names expected by the pipeline):
  - sample_id
  - dapi_path
  - mask_path
  - seg_path
  - reads_unmasked_path
  - reads_masked_path

Example sample sheet content:
```pgsql
sample_id,dapi_path,mask_path,seg_path,reads_unmasked_path,reads_masked_path
E12.5_9,data/dapi/DAPI_E12.5_9.tif,data/masks/DAPI_E12.5_9_mask.tif,data/cellpose/DAPI_E12.5_9_seg.npy,data/unmasked_reads/E12.5_9.csv,data/masked_reads/E12.5_9_masked.csv
```
Checklist before running:
- [ ] All paths listed in the sample sheet exist.
- [ ] Filenames follow the naming convention above.
- [ ] No duplicate sample_id entries.

If you need to rename files to match these patterns, use shell commands to preview then rename (do a dry run first).
