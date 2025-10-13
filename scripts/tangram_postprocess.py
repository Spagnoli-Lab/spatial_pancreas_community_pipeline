#!/usr/bin/env python3
"""Post-Tangram visualisations and summary reports."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc  # type: ignore
import seaborn as sns  # type: ignore
import tangram as tg  # type: ignore
import squidpy as sq  # type: ignore
from matplotlib.colors import ListedColormap


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate plots and summaries from Tangram outputs")
    parser.add_argument("--sample", required=True, help="Sample identifier used for output naming")
    parser.add_argument("--tangram_map", required=True, help="Tangram mapping AnnData (.h5ad) from run_tangram")
    parser.add_argument("--spatial_adata", required=True, help="Tangram-annotated AnnData (.h5ad)")
    parser.add_argument("--output_dir", required=True, help="Directory where visualisations and reports are stored")
    parser.add_argument("--dpi", type=int, default=150, help="Resolution for saved figures")
    return parser.parse_args()


def plot_training_scores(ad_map, output_path: Path, dpi: int) -> None:
    plt.figure()
    tg.plot_training_scores(ad_map, bins=10, alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi)
    plt.close()


def plot_spatial_scatter(adata: sc.AnnData, output_path: Path, dpi: int) -> None:
    cmap = ListedColormap(sns.color_palette("colorblind", 14).as_hex())
    sq.pl.spatial_scatter(
        adata,
        color="Tangram_annotation",
        shape=None,
        figsize=(12, 12),
        size=200,
        spatial_key="spatial",
        cmap=cmap,
        title="Tangram Annotation",
    )
    plt.grid(False)
    handles, labels = plt.gca().get_legend_handles_labels()
    if handles and labels:
        order = list(range(len(labels)))
        if len(order) >= 10:
            order = [4, 5, 6, 7, 8, 9, 0, 1, 2, 3] + [idx for idx in order if idx >= 10]
        plt.legend([handles[i] for i in order], [labels[i] for i in order], bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()


def write_celltype_counts(adata: sc.AnnData, output_path: Path) -> None:
    counts = adata.obs.groupby("Tangram_annotation").size().sort_values(ascending=False)
    counts.to_csv(output_path, sep="\t", header=["count"])


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ad_map = sc.read_h5ad(args.tangram_map)
    spatial_data = sc.read_h5ad(args.spatial_adata)

    plot_training_scores(ad_map, output_dir / f"{args.sample}_tangram_training_scores.png", args.dpi)
    plot_spatial_scatter(spatial_data, output_dir / f"{args.sample}_tangram_scatter.png", args.dpi)
    write_celltype_counts(spatial_data, output_dir / f"{args.sample}_tangram_celltype_counts.txt")

    print(f"[tangram-post] Generated reports under {output_dir}")


if __name__ == "__main__":
    main()
