# Stereo-seq Usage Guide

## Overview  
This guide describes preprocessing and running Polyomino on Stereo-seq spatial transcriptomics data.

## Data Preparation

1. **Obtain processed Stereo-seq data**  
   Use the Stereo-seq pipeline to generate gene count matrices and spatial coordinate files. The data typically includes high-resolution spatial bins.

2. **Load data into AnnData**  
   Make sure the spatial coordinates are in microns and properly formatted in `.obsm['spatial']`. Example:
   ```python
   import scanpy as sc
   adata = sc.read_h5ad('path_to_stereo_seq_data.h5ad')

3. **Quality control**
    Remove low-quality bins if necessary (e.g., bins with fewer than 200 genes detected).

## Generating Spatial Grid

- For Stereo-seq, bin sizes can be very small due to ultra-high resolution.

- Common practice is to aggregate original bins to larger voxel sizes such as **bin5**, **bin25**, or **bin50** to balance resolution and computational load.

- Use the `generate_grid` function with appropriate `width`:

  ```python
  import Polyomino as po
  stdata_grid = po.generate_grid(adata, width=5)  # For bin50 resolution
  ```

- Note: The `width` corresponds to the number of original bins merged in each dimension.

## Running Polyomino

- Initialize Polyomino with the single-cell data and Stereo-seq grid:

  ```python
  po_obj = po.Polyomino(scdata, stdata_grid, cluster_time=1, device='cuda')
  po_obj.allocate()
  ```