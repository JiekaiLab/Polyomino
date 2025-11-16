import argparse

def argument_parser():
    parser = argparse.ArgumentParser(description="CellChIP is an algorithm framework "
                                                 "employing multi-layered regional "
                                                 "constraints to accurately assign "
                                                 "cell locations, enhancing spatial "
                                                 "accuracy and resilience to noise.")

    # Required arguments
    required = parser.add_argument_group('Required')
    required.add_argument("-sc", help="Path of scRNA-Seq anndata", type=str, default=None, required=True)
    required.add_argument("-st", help="Path of spatial anndata", type=str, default=None, required=True)
    required.add_argument("-w", help="width of gridding spatial anndata", type=int, default=None, required=True)

    # Output arguments
    parser.add_argument("-o", "--output", help="Output file", type=str, default="./zmap_assign.h5ad")

    # Custom options
    parser.add_argument("--cluster_time", help="Time allowed for clustering", type=int, default=0)
    parser.add_argument("--custom_region", help="Custom region label for ST data", type=str, default=None)
    parser.add_argument("--cluster_thres", help="Threshold for clustering", type=float, default=None)
    parser.add_argument("--thres", help="Threshold for voxel mapping", type=float, default=0.1)
    parser.add_argument("--method", help="Method used to assing cells to voxels", type=str, choices=['max', 'lap'], default='max')
    parser.add_argument("--device", help="Device to run the computation on", type=str, choices=['cpu', 'cuda'], default='cpu')

    arguments = parser.parse_args()

    return arguments.__dict__
