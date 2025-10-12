#!/usr/bin/env python3
"""Run the pciSeq fit routine using the same inputs prepared in the segmentation notebook.

This utility replicates the notebook step:

    pciSeq.attach_to_log()
    opts = {"save_data": False, "Inefficiency": 0.02}
    res = pciSeq.fit(spots_df, label_image1, scRNA_df2, opts)

By providing CLI arguments, it loads the masked read coordinates, converts the
processed segmentation mask into a sparse matrix, prepares the single-cell
reference profiles, and executes ``pciSeq.fit``. The resulting object can be
pickled for downstream use.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc  # type: ignore
from scipy.sparse import coo_matrix

import pciSeq  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Execute pciSeq.fit on prepared spatial inputs")
    parser.add_argument("--sample", required=True, help="Sample identifier for logging and outputs")
    parser.add_argument("--masked_reads_csv", required=True, help="CSV containing masked read coordinates")
    parser.add_argument(
        "--mask_array",
        required=True,
        help="Path to processed cell mask (.npz with sparse components or dense .npy/.npz array)",
    )
    parser.add_argument(
        "--sc_h5ad",
        required=True,
        help="Precomputed single-cell reference AnnData object (.h5ad)",
    )
    parser.add_argument(
        "--stage",
        help="Developmental stage label to select within the AnnData obs (matches stage_key column)",
    )
    parser.add_argument(
        "--stage_key",
        default="orig.ident",
        help="AnnData obs column name that stores stage identifiers",
    )
    parser.add_argument(
        "--cluster_key",
        default="active.ident",
        help="AnnData obs column containing cluster or cell-type labels",
    )
    parser.add_argument(
        "--exclude_clusters",
        nargs="*",
        default=["Blood cells"],
        help="Cluster labels to remove before running pciSeq (matches values in cluster_key column)",
    )
    parser.add_argument(
        "--collapse_immune",
        action="store_true",
        help="Collapse immune-related cluster labels (T/B/Granulocytes/Macrophages) into 'Immune'",
    )
    parser.add_argument(
        "--inefficiency",
        type=float,
        default=0.02,
        help="Inefficiency parameter to pass inside the opts dictionary",
    )
    parser.add_argument(
        "--save_data",
        action="store_true",
        help="Set opts['save_data']=True so pciSeq writes intermediate data if supported",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where the pickled pciSeq result will be stored",
    )
    parser.add_argument(
        "--result_name",
        default=None,
        help="Optional filename for the pickled result (default: <sample>_pciseq_result.pkl)",
    )
    parser.add_argument(
        "--log_file",
        default=None,
        help="Optional path passed to pciSeq.attach_to_log(); defaults to pciSeq internal behaviour",
    )
    return parser.parse_args()


def load_spots(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={"gene": "Gene", "X": "x", "Y": "y"})
    required = {"Gene", "x", "y"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Masked reads file {path} is missing required columns: {', '.join(sorted(missing))}")
    return df[list(required)]


def load_mask(path: Path) -> coo_matrix:
    if path.suffix == ".npz":
        data = np.load(path, allow_pickle=True)
        if {"data", "row", "col", "shape"}.issubset(data.files):
            return coo_matrix((data["data"], (data["row"], data["col"])), shape=tuple(data["shape"]))
        if "arr_0" in data.files:
            array = data["arr_0"]
            return coo_matrix(array)
        raise ValueError(
            f"Could not interpret NPZ mask at {path}. Expected sparse keys "
            "'data', 'row', 'col', 'shape' or a single array entry."
        )
    if path.suffix == ".npy":
        array = np.load(path, allow_pickle=True)
        return coo_matrix(array)
    raise ValueError(f"Unsupported mask format '{path.suffix}'. Use .npz or .npy files.")


IMMUNE_LABELS = ["T-Cell", "B-Cells", "Granulocytes", "Macrophages"]


def load_sc_reference(
    h5ad_path: Path,
    sample: str,
    stage: str | None,
    stage_key: str,
    cluster_key: str,
    exclude_clusters: list[str],
    collapse_immune: bool,
) -> pd.DataFrame:
    """Return the scRNA reference matrix mimicking scRNA_df2 from the notebook."""
    adata = sc.read_h5ad(h5ad_path)
    stage_value = stage

    if stage_key not in adata.obs.columns:
        if stage_value is not None:
            raise KeyError(f"AnnData.obs is missing stage column '{stage_key}' required for filtering")
        print(f"[pciseq] Stage column '{stage_key}' not found; using all {adata.n_obs} cells.")
    else:
        stage_series = adata.obs[stage_key].astype(str)
        unique_values = set(stage_series.unique())
        if stage_value is None:
            primary = sample.split("_")[0]
            candidates = [primary]
            if "." in primary:
                candidates.append(primary.split(".", 1)[0])
            if primary.replace(".", "") != primary:
                candidates.append(primary.replace(".", ""))
            for candidate in candidates:
                if candidate in unique_values:
                    stage_value = candidate
                    break

        if stage_value is not None:
            mask = stage_series == stage_value
            if not mask.any():
                available = ", ".join(sorted(unique_values))
                raise ValueError(
                    f"No cells found for stage '{stage_value}' in column '{stage_key}'. "
                    f"Available values: {available}"
                )
            adata = adata[mask].copy()
            print(f"[pciseq] Selected {mask.sum()} cells for stage '{stage_value}' (column '{stage_key}').")
        else:
            available = ", ".join(sorted(unique_values))
            print(
                "[pciseq] No matching stage found for sample "
                f"'{sample}' (available stages: {available}); using all {adata.n_obs} cells."
            )

    if cluster_key not in adata.obs.columns:
        raise KeyError(f"AnnData.obs is missing cluster column '{cluster_key}' required for pciSeq input")

    if exclude_clusters:
        cluster_series = adata.obs[cluster_key].astype(str)
        mask = ~cluster_series.isin(exclude_clusters)
        adata = adata[mask].copy()
        cluster_series = adata.obs[cluster_key].astype(str)
        dropped = (~mask).sum()
        print(f"[pciseq] Excluded clusters {exclude_clusters}; removed {dropped} cells, remaining {adata.n_obs}.")
    else:
        cluster_series = adata.obs[cluster_key].astype(str)

    if collapse_immune:
        collapsed = cluster_series.replace(IMMUNE_LABELS, "Immune")
        adata.obs[cluster_key] = collapsed
        cluster_series = collapsed
        print("[pciseq] Collapsed immune-related clusters into 'Immune'.")

    if adata.n_obs == 0:
        raise ValueError("No cells remain in the single-cell reference after filtering.")

    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    reference_df = pd.DataFrame(matrix.T, index=adata.var_names, columns=adata.obs_names)
    rename_map = cluster_series.to_dict()
    reference_df.rename(columns=rename_map, inplace=True)
    return reference_df


def main() -> None:
    args = parse_args()

    sample = args.sample
    masked_reads_path = Path(args.masked_reads_csv)
    mask_path = Path(args.mask_array)
    sc_h5ad_path = Path(args.sc_h5ad)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    spots_df = load_spots(masked_reads_path)
    label_sparse = load_mask(mask_path)
    sc_reference = load_sc_reference(
        sc_h5ad_path,
        sample=sample,
        stage=args.stage,
        stage_key=args.stage_key,
        cluster_key=args.cluster_key,
        exclude_clusters=args.exclude_clusters,
        collapse_immune=args.collapse_immune,
    )

    if args.log_file:
        pciSeq.attach_to_log(args.log_file)
    else:
        pciSeq.attach_to_log()

    opts = {"save_data": args.save_data, "Inefficiency": args.inefficiency}
    res = pciSeq.fit(spots_df, label_sparse, scRNAseq=sc_reference, opts=opts)

    result_name = args.result_name or f"{sample}_pciseq_result.pkl"
    result_path = output_dir / result_name
    with open(result_path, "wb") as handle:
        pickle.dump(res, handle)

    print(f"[pciseq] Completed fit for sample '{sample}'. Result saved to: {result_path}")


if __name__ == "__main__":
    main()
