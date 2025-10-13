# Spatial Pipeline Flow (Beginner Friendly)

```mermaid
flowchart TD
    A[samplesheet.csv] --> B[Normalize sample rows<br/>(channel enrichment)]
    B --> C[FILTER_CELLS_BY_TISSUE<br/>cell_segmentation_filter.py]
    C -->|processed mask + viz + masked reads path| D[RUN_PCISEQ<br/>run_pciseq_fit.py]
    D -->|filtered tuple| E[OVERLAY_MASKED_READS<br/>overlay_masked_reads.py]
    E -->|"QC overlays"<br/>overlay_png_summary]

    D -->|pciSeq result pickle| F[PCISEQ_TO_TANGRAM_INPUT<br/>pciseq_to_anndata.py]
    F -->|pciSeq AnnData| G[RUN_TANGRAM<br/>run_tangram.py]
    G -->|spatial AnnData + projected genes + map| H[TANGRAM_PLOT<br/>tangram_postprocess.py]

    classDef data fill:#f3f4ff,stroke:#555,color:#222;
    classDef process fill:#e0f7f4,stroke:#20706f,color:#0f3d39;
    class A data;
    class B,C,D,E,F,G,H process;
```

**Legend**

- **samplesheet.csv** – provides per-sample paths for masks, segmentations, and masked reads.
- **FILTER_CELLS_BY_TISSUE** – rebuilds/filters the segmentation mask.
- **RUN_PCISEQ** – overlays pciSeq classifications on the filtered mask.
- **OVERLAY_MASKED_READS** – generates QC figures with masked reads on the mask.
- **PCISEQ_TO_TANGRAM_INPUT** – converts pciSeq output into an AnnData compatible with Tangram.
- **RUN_TANGRAM** – aligns single-cell reference to spatial data, producing annotated AnnData plus projected genes and a mapping object.
- **TANGRAM_PLOT** – produces training-score plots, annotated scatter, and a cell-type count report.

Outputs are written under `outputs/nf/` in folders matching each major stage (e.g. `cell_filter`, `overlay_qc`, `pciseq`, `tangram`, `tangram/post`).
