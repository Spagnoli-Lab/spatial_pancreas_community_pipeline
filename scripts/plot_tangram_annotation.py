#!/usr/bin/env python3
"""Generate Tangram cluster annotation plots from spatial AnnData outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc  # type: ignore
import tangram as tg  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Tangram cluster annotations for a spatial sample")
    parser.add_argument("--sample", required=True, help="Sample identifier used for naming outputs")
    parser.add_argument("--sc_reference", required=True, help="Path to single-cell reference AnnData (.h5ad)")
    parser.add_argument("--tangram_anndata", required=True, help="Tangram-annotated AnnData (.h5ad)")
    parser.add_argument("--cluster_key", default="active.ident", help="Column in sc_reference.obs with cluster labels")
    parser.add_argument("--output_dir", required=True, help="Directory where the figure will be saved")
    parser.add_argument("--spot_size", type=float, default=45.0, help="Scatter size for spatial plot")
    parser.add_argument("--dpi", type=int, default=150, help="Resolution for the saved figure")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    scdata = sc.read_h5ad(args.sc_reference)
    spatial_data = sc.read_h5ad(args.tangram_anndata)

    if args.cluster_key not in scdata.obs.columns:
        available = ", ".join(scdata.obs.columns)
        raise KeyError(
            f"Cluster column '{args.cluster_key}' missing from sc reference AnnData.obs. "
            f"Available columns: {available}"
        )

    annotation_list = list(pd.unique(scdata.obs[args.cluster_key]))

    tg.plot_cell_annotation_sc(
        spatial_data,
        annotation_list,
        perc=0.02,
        scale_factor=1,
        spot_size=args.spot_size,
        x="x",
        y="y",
    )

    fig = plt.gcf()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    figure_path = output_dir / f"{args.sample}_tangram_annotation.png"
    fig.savefig(figure_path, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    print(f"[tangram-plot] Saved Tangram cluster plot to {figure_path}")


if __name__ == "__main__":
    main()
