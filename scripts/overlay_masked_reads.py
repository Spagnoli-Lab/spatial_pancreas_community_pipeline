#!/usr/bin/env python3
"""Overlay masked read coordinates onto processed cell segmentation masks.

This utility is intended to be run within the Nextflow pipeline after
cell segmentation filtering. It reconstructs the processed mask from the
compressed sparse representation, converts it to a color label image, and
plots the masked read coordinates on top as a QC visualization.

Usage example:
    python scripts/overlay_masked_reads.py \\
        --sample E14.5_5 \\
        --masked_reads_csv data/masked_reads/E14.5_5_reads.csv \\
        --mask_npz data/masks/E14.5_5_mask.npz \\
        --output_dir results/overlays

You may instead provide ``--sample_stage`` and ``--mouse_index`` if that fits
your pipeline inputs better. The script writes ``*_masked_reads_overlay.png``
and a tab-delimited ``*_masked_reads_summary.txt`` into ``--output_dir``.
Pass ``--show`` while running locally to pop up the overlay window for a quick
visual QC that masks and coordinates align.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from skimage.color import label2rgb


def parse_sample(sample_id: str) -> tuple[str, str]:
    if '_' not in sample_id:
        raise ValueError(
            f"Sample identifier '{sample_id}' must include an underscore separating stage and index"
        )
    stage_part, index_part = sample_id.rsplit('_', 1)
    if not stage_part or not index_part:
        raise ValueError(
            f"Sample identifier '{sample_id}' must provide both stage and index"
        )
    return stage_part, index_part


def resolve_sample_identifier(sample: str | None, stage: str | None, index: str | None) -> tuple[str, str, str]:
    """Return normalized sample identifier and its stage/index components.

    Users may provide the full sample id (e.g. 'E14.5_5') or split pieces.
    The function ensures the pieces are consistent and always returns the
    resolved sample id alongside stage and index strings.
    """
    stage = None if stage is None else str(stage)
    index = None if index is None else str(index)

    if sample:
        sample_id = sample
        stage_part, index_part = parse_sample(sample_id)
        if stage and stage_part != stage:
            raise ValueError(
                f"Provided stage '{stage}' does not match sample id stage '{stage_part}'"
            )
        if index and index_part != index:
            raise ValueError(
                f"Provided mouse index '{index}' does not match sample id index '{index_part}'"
            )
    else:
        if not stage or not index:
            raise ValueError(
                "Provide either --sample or both --sample_stage and --mouse_index"
            )
        sample_id = f"{stage}_{index}"
        stage_part, index_part = stage, index

    return sample_id, stage_part, index_part


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Overlay masked read coordinates onto cell masks")
    parser.add_argument("--sample", help="Full sample identifier, e.g. 'E14.5_5'")
    parser.add_argument("--sample_stage", help="Sample developmental stage, e.g. 'E14.5'")
    parser.add_argument("--mouse_index", help="Mouse index within a stage, e.g. '5'")
    parser.add_argument("--masked_reads_csv", required=True, help="Path to masked reads CSV for the sample")
    parser.add_argument("--mask_npz", required=True, help="Path to processed mask stored as a compressed sparse matrix (.npz)")
    parser.add_argument("--output_dir", required=True, help="Directory where outputs will be written")
    parser.add_argument("--dpi", type=int, default=72, help="Resolution used for the saved figure")
    parser.add_argument("--point_size", type=float, default=10.0, help="Scatter point size for masked reads")
    parser.add_argument("--point_color", default="white", help="Scatter point color")
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the overlay interactively after saving for quick visual QC",
    )
    return parser.parse_args()


def load_masked_reads(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=",")
    # Standardize expected column names
    df = df.rename(columns={"gene": "Gene", "X": "x", "Y": "y"})

    required_columns = {"Gene", "x", "y"}
    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(f"Masked reads file {path} is missing required columns: {', '.join(sorted(missing))}")

    return df[["Gene", "x", "y"]]


def load_mask(npz_path: Path) -> np.ndarray:
    data = np.load(npz_path)
    if not {"data", "row", "col", "shape"}.issubset(data.files):
        raise ValueError(f"Processed mask file {npz_path} is missing expected sparse matrix keys")

    mask_sparse = coo_matrix((data["data"], (data["row"], data["col"])), shape=data["shape"])
    return mask_sparse.toarray()


def create_overlay(
    label_image: np.ndarray,
    spots_df: pd.DataFrame,
    sample: str,
    output_path: Path,
    dpi: int,
    point_size: float,
    point_color: str,
    show: bool,
) -> None:
    rgb_label_image = label2rgb(label_image, bg_label=0)

    height, width = label_image.shape
    fig_width = width / dpi
    fig_height = height / dpi

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=dpi)
    ax.imshow(rgb_label_image)
    ax.scatter(spots_df["x"], spots_df["y"], s=point_size, c=point_color)
    ax.axis("off")
    fig.tight_layout(pad=0)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    masked_reads_path = Path(args.masked_reads_csv)
    mask_npz_path = Path(args.mask_npz)

    sample_id, stage_part, index_part = resolve_sample_identifier(
        args.sample, args.sample_stage, args.mouse_index
    )

    print(f"[overlay] Sample: {sample_id}")
    print(f"[overlay] Stage: {stage_part}")
    print(f"[overlay] Mouse index: {index_part}")

    print(f"[overlay] Reading masked reads from: {masked_reads_path}")
    spots_df = load_masked_reads(masked_reads_path)
    print(f"[overlay] Loaded {len(spots_df)} masked reads")

    print(f"[overlay] Loading processed mask from: {mask_npz_path}")
    label_image = load_mask(mask_npz_path)
    print(f"[overlay] Mask dimensions: {label_image.shape}")

    overlay_path = output_dir / f"{sample_id}_masked_reads_overlay.png"
    print(f"[overlay] Saving overlay figure to: {overlay_path}")

    create_overlay(
        label_image=label_image,
        spots_df=spots_df,
        sample=sample_id,
        output_path=overlay_path,
        dpi=args.dpi,
        point_size=args.point_size,
        point_color=args.point_color,
        show=args.show,
    )

    summary_path = output_dir / f"{sample_id}_masked_reads_summary.txt"
    print(f"[overlay] Writing summary to: {summary_path}")
    with open(summary_path, "w", encoding="utf-8") as fh:
        fh.write(f"sample\tstage\tmouse_index\tmasked_reads\theight\twidth\n")
        fh.write(
            f"{sample_id}\t{stage_part}\t{index_part}\t{len(spots_df)}\t{label_image.shape[0]}\t{label_image.shape[1]}\n"
        )

    print("[overlay] Completed overlay generation")


if __name__ == "__main__":
    main()
