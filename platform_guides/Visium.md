# Visium Usage Guide

## Overview
This guide explains how to preprocess 10x Genomics Visium spatial transcriptomics data for use with Polyomino.

## Data Preparation

1. **Obtain processed Visium data**  
   Use Space Ranger to generate filtered feature-barcode matrices (`filtered_feature_bc_matrix.h5`) and spatial coordinate files.

2. **Load data into AnnData**  
   Use Scanpy or other tools to read Visium data:
   ```python
   import scanpy as sc
   adata = sc.read_visium('path_to_visium_folder')

   3.**Quality control**
 	Filter spots or genes as needed to remove low-quality data.

## Generating Spatial Grid

- Polyomino uses voxelization to create spatial bins:

  - Typical voxel size for Visium: 1
  - Adjust `width` depending on tissue resolution and analysis goals

- Generate the grid:

  ```python
  import Polyomino as po
  stdata_grid = po.generate_grid(adata, width=1)
  ```

## Running Polyomino

- Initialize and run mapping:

  ```python
  po_obj = po.Polyomino3D(scdata, stdata_grid, cluster_time=1, device='cpu')
  po_obj.allocate()
  ```

