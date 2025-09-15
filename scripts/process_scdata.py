import argparse
import pandas as pd
import scanpy as sc

def process_scdata():
    parser = argparse.ArgumentParser(description="Preprocess single-cell data for spatial analysis.")
    parser.add_argument('--scRNA_df', required=True, help='Path to single-cell expression matrix (TXT)')
    parser.add_argument('--scRNA_meta_df', required=True, help='Path to single-cell metadata (TXT)')
    parser.add_argument('--stage', default='E14', help='Developmental stage to filter (default: E14)')
    parser.add_argument('--output_h5ad', required=True, help='Output path and file name for processed AnnData object (add .h5ad to the end)')
    args = parser.parse_args()

    print("=== Starting Single-Cell Data Preprocessing ===")

    # Load data
    print("Loading single-cell expression data...")
    scRNA_df = pd.read_csv(args.scRNA_df, sep=" ", index_col=0)
    print(f"Loaded expression data: {scRNA_df.shape[0]} genes x {scRNA_df.shape[1]} cells")

    print("Loading single-cell metadata...")
    scRNA_meta_df = pd.read_csv(args.scRNA_meta_df, sep=" ", index_col=0)
    print(f"Loaded metadata: {scRNA_meta_df.shape[0]} cells x {scRNA_meta_df.shape[1]} attributes")

    # Transpose the single cell data to have genes as columns
    print("Transposing expression matrix...")
    scRNA_df = scRNA_df.transpose()

    # Change the stage
    print(f"Filtering cells by developmental stage: {args.stage}")
    scRNA_df = scRNA_df[scRNA_meta_df["orig.ident"] == args.stage]
    print(f"After stage filtering: {scRNA_df.shape[0]} cells")

    # Filter out blood cells
    print("Filtering out blood cells...")
    scRNA_df = scRNA_df[(scRNA_meta_df["active.ident"] != "Blood cells")]
    scRNA_meta_df = scRNA_meta_df[scRNA_meta_df["active.ident"] != "Blood cells"]
    print(f"After blood cell filtering: {scRNA_df.shape[0]} cells")

    # Replace T-Cell, B-Cells, Granulocytes, Macrophages with Immune
    print("Merging immune cell types...")
    scRNA_meta_df.replace(["T-Cell", "B-Cells", "Granulocytes", "Macrophages"], "Immune", inplace=True)
    print("Immune cell types merged into 'Immune' category")

    # Transpose the single cell data to have cells as columns
    scRNA_df = scRNA_df.transpose()
    
    # Rename the columns to the cell type labels
    print("Renaming columns to cell type labels...")
    scRNA_df = scRNA_df.rename(columns=scRNA_meta_df["active.ident"])

    # Create single-cell AnnData object
    print("Creating AnnData object...")
    scdata = sc.AnnData(scRNA_df.transpose())  # Transpose so genes are columns
    scdata.obs["cluster_label"] = scRNA_df.columns
    scdata.var.index = scRNA_df.index
    sc.pp.normalize_total(scdata)
    print(f"Created AnnData object: {scdata.shape[0]} cells x {scdata.shape[1]} genes")

    # Save the processed AnnData object
    print(f"Saving processed data to: {args.output_h5ad}")
    scdata.write(args.output_h5ad)
    print("=== Single-Cell Data Preprocessing Complete ===")

if __name__ == "__main__":
    process_scdata()