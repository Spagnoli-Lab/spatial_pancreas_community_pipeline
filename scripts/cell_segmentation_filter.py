#!/usr/bin/env python3
"""
Cell Segmentation and Tissue Mask Filtering Pipeline

This script processes cellpose segmentation results and applies tissue mask filtering
to remove cells that fall outside the tissue boundary. It loads cellpose segmentation
files, applies tissue mask filtering, and generates both processed masks and visualizations.

Usage:
    python cell_segmentation_filter.py --sample E14.5_9 --cellpose_path /path/to/cellpose --mask_path /path/to/masks --output_dir /path/to/output

Inputs:
    - Cellpose segmentation file (numpy .npy format)
    - Tissue mask image (.tif format)
    - Sample name for file identification

Outputs:
    - Processed cell mask (sparse matrix format)
    - Visualization of filtered cell masks (.tif format)

Author: Siwanart Ma
Date: 13/09/2025
"""

import argparse
import os
import numpy as np
from PIL import Image
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
import skimage.color
from skimage import io as skio
from skimage.measure import regionprops
import re

def cell_segmentation_filter():
    """
    Main function to process cellpose segmentation and apply tissue mask filtering.
    
    This function:
    1. Loads cellpose segmentation results
    2. Loads tissue mask image
    3. Filters out cells whose centroids fall outside the tissue mask
    4. Saves processed mask as sparse matrix
    5. Generates visualization of filtered cells
    """
    parser = argparse.ArgumentParser(description="Process cellpose masks and apply tissue mask filtering.")
    parser.add_argument('--sample', required=True, help='Sample name (e.g., E14.5_9)')
    parser.add_argument('--cellpose_path', required=True, help='Path to cellpose files directory')
    parser.add_argument('--mask_path', required=True, help='Path to tissue mask files directory')
    parser.add_argument('--output_dir', required=True, help='Output directory for processed masks')
    args = parser.parse_args()

    print("=== Starting Cell Segmentation and Tissue Mask Filtering ===")
    print(f"Processing sample: {args.sample}")

    # Find the cellpose file for this sample
    print("Loading cellpose segmentation file...")
    cellpose_files = [filename for filename in os.listdir(args.cellpose_path) 
                     if filename.startswith(f"DAPI_{args.sample}_")]
    
    if not cellpose_files:
        raise FileNotFoundError(f"No cellpose file found for sample {args.sample}")
    
    cellpose_file = cellpose_files[0]
    print(f"Found cellpose file: {cellpose_file}")

    # Load cellpose data
    dat = np.load(os.path.join(args.cellpose_path, cellpose_file), allow_pickle=True).item()
    label_image = dat['masks']
    print(f"Loaded label image with shape: {label_image.shape}")
    
    # Find mask file using regex pattern
    # Accept either DAPI_<sample>_mask.tif or DAPI_<sample>_<anything>_mask.tif
    exact_pattern = re.compile(f"^DAPI_{re.escape(args.sample)}_mask\\.tif$")
    extended_pattern = re.compile(f"^DAPI_{re.escape(args.sample)}_.*_mask\\.tif$")
    mask_files = [
        f for f in os.listdir(args.mask_path)
        if exact_pattern.match(f) or extended_pattern.match(f)
    ]

    if not mask_files:
        raise FileNotFoundError(
            f"No mask file found for sample {args.sample}. Expected 'DAPI_{args.sample}_mask.tif' or 'DAPI_{args.sample}_*_mask.tif'"
        )

    mask_filename = mask_files[0]
    mask_path = os.path.join(args.mask_path, mask_filename)
    print(f"Found mask file: {mask_filename}")

    # Use skimage to read TIFF (more robust for ImageJ/BigTIFF)
    maskarray = skio.imread(mask_path)
    print(f"Loaded tissue mask with shape: {maskarray.shape}")

    # Process cell regions
    print("Processing cell regions and applying tissue mask filtering...")
    props = regionprops(label_image)
    cells_removed = 0

    for z in props:
        x_pos = int(z['centroid'][0])
        y_pos = int(z['centroid'][1])
        coords = z["coords"]
        
        # Check if cell centroid is outside tissue mask
        if maskarray[x_pos, y_pos] == 0:
            # Remove this cell by setting all its pixels to 0
            for n in coords:
                label_image[n[0]][n[1]] = 0
            cells_removed += 1

    print(f"Removed {cells_removed} cells outside tissue mask")

    # Convert to sparse matrix
    print("Converting to sparse matrix format...")
    label_image1 = coo_matrix(label_image)

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Save processed mask as sparse matrix
    mask_output_path = os.path.join(args.output_dir, f"{args.sample}_processed_mask.npz")
    np.savez_compressed(mask_output_path, 
                       data=label_image1.data,
                       row=label_image1.row,
                       col=label_image1.col,
                       shape=label_image1.shape)
    print(f"Saved processed mask to: {mask_output_path}")

    # Generate visualization
    print("Generating visualization...")
    rgb_label_image = skimage.color.label2rgb(label_image1.toarray(), bg_label=0)
    
    # Set up matplotlib for saving
    _dpi = 72
    plt.figure(figsize=(1500/_dpi, 1500/_dpi), dpi=_dpi)
    imgplot = plt.imshow(rgb_label_image)
    plt.grid(False)
    
    # Save visualization
    viz_output_path = os.path.join(args.output_dir, f"{args.sample}_cell_masks.tif")
    plt.savefig(viz_output_path, dpi=_dpi, bbox_inches='tight')
    plt.close()  # Close the figure to free memory
    print(f"Saved visualization to: {viz_output_path}")

    print("=== Cell Segmentation and Tissue Mask Filtering Complete ===")

if __name__ == "__main__":
    cell_segmentation_filter()
