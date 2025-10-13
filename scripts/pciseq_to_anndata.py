#!/usr/bin/env python3
"""Convert pciSeq fit results into an AnnData object ready for Tangram."""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scanpy as sc  # type: ignore


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Transform pciSeq result pickle into an AnnData .h5ad file")
    parser.add_argument("--sample", required=True, help="Sample identifier used for output naming")
    parser.add_argument(
        "--pciseq_result",
        required=True,
        help="Path to the pciSeq result pickle produced by run_pciseq_fit.py",
    )
    parser.add_argument("--output_dir", required=True, help="Directory where the AnnData file will be written")
    parser.add_argument(
        "--immune_labels",
        nargs="*",
        default=["T-Cell", "B-Cells", "Granulocytes", "Macrophages"],
        help="Cluster labels already collapsed into 'Immune'; used only for logging",
    )
    return parser.parse_args()


def _ensure_iterable(value: Iterable) -> list:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return [value]


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.pciseq_result, "rb") as handle:
        res = pickle.load(handle)

    if not isinstance(res, (list, tuple)) or len(res) < 2:
        raise ValueError(f"Unexpected pciSeq result structure in {args.pciseq_result!r}")

    cells_df = pd.DataFrame(res[0])
    genes_df = pd.DataFrame(res[1])

    if "Gene" not in genes_df.columns:
        raise KeyError("pciSeq genes dataframe is missing 'Gene' column")

    genes = pd.Index(pd.unique(genes_df["Gene"]))
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}

    n_cells = len(cells_df)
    n_genes = len(genes)
    counts = np.zeros((n_cells, n_genes), dtype=float)

    genenames_series = cells_df.get("Genenames")
    cellgene_series = cells_df.get("CellGeneCount")
    if genenames_series is None or cellgene_series is None:
        raise KeyError("pciSeq cells dataframe must include 'Genenames' and 'CellGeneCount' columns")

    for row_idx, (gene_names, gene_counts) in enumerate(zip(genenames_series, cellgene_series)):
        for gene, count in zip(_ensure_iterable(gene_names), _ensure_iterable(gene_counts)):
            if gene in gene_to_idx:
                counts[row_idx, gene_to_idx[gene]] = count

    class_series = cells_df.get("ClassName")
    prob_series = cells_df.get("Prob")
    celltypes: list[str] = []
    if class_series is not None and prob_series is not None:
        for class_names, probs in zip(class_series, prob_series):
            names = _ensure_iterable(class_names)
            probabilities = _ensure_iterable(probs)
            if names and probabilities:
                max_idx = int(np.argmax(probabilities))
                celltypes.append(str(names[max_idx]))
            elif names:
                celltypes.append(str(names[0]))
            else:
                celltypes.append("Unknown")
    elif class_series is not None:
        celltypes = [str(_ensure_iterable(names)[0]) if _ensure_iterable(names) else "Unknown" for names in class_series]
    else:
        celltypes = ["Unknown"] * n_cells

    spatial_cols = [col for col in ("X", "Y") if col in cells_df.columns]
    if len(spatial_cols) != 2:
        raise KeyError("pciSeq cells dataframe must include 'X' and 'Y' columns for spatial coordinates")
    spatial = cells_df[["X", "Y"]].to_numpy(dtype=float, copy=True)

    obs_index = pd.Index([str(i) for i in range(n_cells)])
    obs = pd.DataFrame(
        {
            "celltype": celltypes,
            "x": spatial[:, 0],
            "y": spatial[:, 1],
        },
        index=obs_index,
    )

    var = pd.DataFrame(index=genes, data={"gene_ids": genes})

    adata = sc.AnnData(X=counts, obs=obs, var=var)
    adata.var_names_make_unique()
    adata.obsm["spatial"] = spatial

    sc.pp.normalize_total(adata)

    output_path = output_dir / f"{args.sample}_pciseq_predicted.h5ad"
    adata.write(output_path)

    print(f"[pciseq->anndata] Wrote AnnData file to {output_path}")


if __name__ == "__main__":
    main()
