#!/usr/bin/env python3
"""Run Tangram to map scRNA reference labels onto pciSeq-derived spatial data."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc  # type: ignore
import tangram as tg  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Execute Tangram using preprocessed AnnData inputs")
    parser.add_argument("--sample", required=True, help="Sample identifier used for output naming")
    parser.add_argument("--sc_reference", required=True, help="Path to single-cell AnnData reference (.h5ad)")
    parser.add_argument("--pciseq_anndata", required=True, help="Path to AnnData generated from pciSeq results")
    parser.add_argument("--cluster_key", default="active.ident", help="Column in sc_reference.obs with cluster labels")
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where Tangram artifacts (h5ad + figures) will be written",
    )
    parser.add_argument(
        "--num_epochs",
        type=int,
        default=1000,
        help="Number of training epochs for Tangram mapping",
    )
    parser.add_argument("--device", default="cpu", help="Compute device for Tangram (e.g. 'cpu' or 'cuda:0')")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    scdata = sc.read_h5ad(args.sc_reference)
    spatial_data = sc.read_h5ad(args.pciseq_anndata)

    spatial_data.obs_names = spatial_data.obs_names.astype(str)
    scdata.obs_names = scdata.obs_names.astype(str)

    genes = scdata.var.index.to_list()
    tg.pp_adatas(scdata, spatial_data, genes=genes)

    ad_map = tg.map_cells_to_space(
        scdata,
        spatial_data,
        device=args.device,
        num_epochs=args.num_epochs,
        mode="cells",
        cluster_label=args.cluster_key,
        scale=True
    )

    ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=scdata)

    tg.project_cell_annotations(ad_map, spatial_data, annotation=args.cluster_key)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    result = spatial_data.obsm["tangram_ct_pred"].idxmax(axis=1)
    spatial_data.obs["Tangram_annotation"] = result.to_numpy()

    tangram_path = output_dir / f"{args.sample}_tangram_predicted.h5ad"
    spatial_data.write(tangram_path)
    print(f"[tangram] Wrote Tangram-annotated AnnData to {tangram_path}")

    projected_path = output_dir / f"{args.sample}_projected.h5ad"
    ad_ge.write(projected_path)
    print(f"[tangram] Wrote projected gene expression to {projected_path}")

    map_path = output_dir / f"{args.sample}_tangram_map.h5ad"
    ad_map.write(map_path)
    print(f"[tangram] Wrote Tangram mapping AnnData to {map_path}")


if __name__ == "__main__":
    main()
