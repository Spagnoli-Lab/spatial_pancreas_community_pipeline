# Python Script Reference

Overview of the helper utilities that back the Nextflow workflow. Each entry explains what the script contributes, the main arguments you can pass on the command line, and the expected inputs and outputs.

## apply_masked_dapi.py
- **What & why**: Rebuilds masked DAPI TIFFs by combining 8-bit tissue masks with their raw fluorescence images. Supports quick QC outside the pipeline and generates both full-frame and cropped variants for downstream review.
- **Key arguments**: `--mask-dir`, `--dapi-dir`, `--output-dir`, `--overwrite`, `--dry-run`, `--verbose`.
- **Inputs**: One TIFF mask per sample plus the matching raw DAPI TIFF.
- **Outputs**: Two TIFFs per sample (`full/…_masked_dapi.tif` and `cut_black_edge/…_masked_dapi.tif`) in the target directory.

## cell_segmentation_filter.py
- **What & why**: Loads Cellpose segmentation results, removes cells falling outside a tissue mask, and emits filtered masks for later steps. Keeps only tissue-supported cells to reduce downstream noise.
- **Key arguments**: `--sample`, `--cellpose_path`, `--mask_path`, `--output_dir`.
- **Inputs**: Cellpose `.npy` segmentation and a tissue mask TIFF for the same sample.
- **Outputs**: Compressed sparse mask (`*_processed_mask.npz`) and a TIFF visualization (`*_cell_masks.tif`).

## overlay_masked_reads.py
- **What & why**: Recreates the filtered mask, overlays masked read coordinates, and writes a QC PNG plus a text summary. Used to verify that filtered segmentation aligns with captured reads.
- **Key arguments**: `--masked_reads_csv`, `--mask_npz`, `--output_dir`, optionally `--sample` or `--sample_stage`/`--mouse_index`, plus plotting controls (`--dpi`, `--point_size`, `--point_color`, `--show`).
- **Inputs**: Processed mask (`*_processed_mask.npz`) and masked reads CSV with `Gene`, `x`, `y` columns.
- **Outputs**: Overlay image (`*_masked_reads_overlay.png`) and tab-delimited summary (`*_masked_reads_summary.txt`).

## run_pciseq_fit.py
- **What & why**: Runs `pciSeq.fit` to assign cell types using masked reads, segmentation, and single-cell reference profiles. Produces the pickled object that downstream steps convert to AnnData.
- **Key arguments**: `--sample`, `--masked_reads_csv`, `--mask_array`, `--sc_h5ad`, `--output_dir`, stage/cluster controls (`--stage`, `--stage_key`, `--cluster_key`, `--exclude_clusters`, `--collapse_immune`), plus run-time options (`--inefficiency`, `--save_data`, `--result_name`, `--log_file`).
- **Inputs**: Masked reads CSV, processed mask (`.npz` or `.npy`), and single-cell reference `.h5ad`.
- **Outputs**: Pickled pciSeq result (`*_pciseq_result.pkl`) in the chosen directory.

## pciseq_to_anndata.py
- **What & why**: Converts the pciSeq pickle into an AnnData object that Tangram expects, normalising counts and storing inferred cell types and spatial coordinates.
- **Key arguments**: `--sample`, `--pciseq_result`, `--output_dir`, optional `--immune_labels` list for reporting.
- **Inputs**: Pickled pciSeq fit result from `run_pciseq_fit.py`.
- **Outputs**: AnnData file (`*_pciseq_predicted.h5ad`) ready for Tangram.

## run_tangram.py
- **What & why**: Executes Tangram to map single-cell labels onto spatial data, producing the annotated spatial AnnData plus the projected gene expression and mapping objects.
- **Key arguments**: `--sample`, `--sc_reference`, `--pciseq_anndata`, `--cluster_key`, `--output_dir`, `--num_epochs`, `--device`.
- **Inputs**: Single-cell reference AnnData and the pciSeq-derived AnnData.
- **Outputs**: Three AnnData artifacts (`*_tangram_predicted.h5ad`, `*_projected.h5ad`, `*_tangram_map.h5ad`) saved under the output directory.

## tangram_postprocess.py
- **What & why**: Generates visual QA (training curves, scatter plot) and a tabular cell-type summary from Tangram outputs. Summaries help flag mapping issues early.
- **Key arguments**: `--sample`, `--tangram_map`, `--spatial_adata`, `--output_dir`, `--dpi`.
- **Inputs**: Tangram mapping `.h5ad` and the annotated spatial `.h5ad` produced by `run_tangram.py`.
- **Outputs**: Training curve PNG, spatial scatter PNG, and cell-type counts TSV written to the output directory.

## plot_tangram_annotation.py
- **What & why**: Optional helper that re-plots Tangram annotations using a custom cluster list from the single-cell reference. Useful for publication-style figures or quick inspection.
- **Key arguments**: `--sample`, `--sc_reference`, `--tangram_anndata`, `--cluster_key`, `--output_dir`, `--spot_size`, `--dpi`.
- **Inputs**: Single-cell reference `.h5ad` and Tangram-annotated `.h5ad`.
- **Outputs**: PNG figure (`*_tangram_annotation.png`) saved in the output directory.

## process_scdata.py
- **What & why**: Standalone utility for preprocessing tabular single-cell data into an AnnData reference by filtering stages, removing blood cells, merging immune labels, and normalising counts. Useful when regenerating the reference input for the pipeline.
- **Key arguments**: `--scRNA_df`, `--scRNA_meta_df`, `--stage`, `--output_h5ad`.
- **Inputs**: Expression matrix TXT and metadata TXT exported from upstream single-cell workflows.
- **Outputs**: Processed single-cell AnnData reference written to the specified `.h5ad`.
