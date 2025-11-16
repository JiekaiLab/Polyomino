from .utilities.util import *
from .model.model import Cell2Clusters, Cell2Spots, Cell2VisiumSpots, Cell2SpotsStrips,Cell2SpotsEnhance
from .model.ot_model import solve_OT
from .Polyomino_argument import argument_parser
import pandas as pd
from scipy.sparse import issparse
import numpy as np
import scanpy as sc
import tqdm
import torch
from sklearn.metrics.pairwise import cosine_similarity as cos_s

class Polyomino:
    """
    This class runs the Polyomino algorithm on scRNA-seq data and spatial transcriptomics or bulk RNA datasets.
    It reconstructs spatial transcriptomics data based on the Polyomino algorithm.

    Attributes:
        scdata (anndata.AnnData): scRNA-seq data.
        stdata (anndata.AnnData): Spatial transcriptomics data.
        bulkX (anndata.AnnData): Bulk transcriptomics data (X-axis).
        bulkY (anndata.AnnData): Bulk transcriptomics data (Y-axis).
        genes (list): List of genes to use for mapping.
        histology (np.array): Histology data.
        cluster_time (int): Number of clustering iterations.
        custom_label (str): Custom label for clustering.
        pca (str): PCA method to use.
        n_pcs (int): Number of principal components.
        pc_frac (float): Fraction of PCs to use.
        target_num (int): Number of clusters.
        cluster_thres (float): Threshold to filter small values in clustering.
        ot_cluster_thres (float): Threshold for OT clustering.
        emptygrid (np.array): Empty grid matrix.
        shape (str): Shape of spatial transcriptomics data ('square' or 'hexagon').
        device (str): Device for training ('cpu' or 'cuda').
        cluster_matrix (np.ndarray): Matrix for clustering results.
        spot_matrix (np.ndarray): Matrix for spot mapping results.
    """

    def __init__(
        self, 
        scdata, 
        stdata=None,
        bulkX=None, 
        bulkY=None, 
        genes=None, 
        histology=None, 
        cluster_time=0, 
        custom_label=None,
        pca=None, 
        n_pcs=50, 
        pc_frac=1.0, 
        target_num=10,
        cluster_thres=None,
        ot_cluster_thres=0.01,
        emptygrid=None,
        shape="square", 
        device='cpu'
    ):
        """
        Initializes the Polyomino class with provided parameters.

        Args:
            scdata (anndata.AnnData): scRNA-seq data.
            stdata (anndata.AnnData): Spatial transcriptomics data.
            bulkX (anndata.AnnData): Bulk transcriptomics data (X-axis).
            bulkY (anndata.AnnData): Bulk transcriptomics data (Y-axis).
            genes (list): List of genes to use for mapping.
            histology (np.array): Histology data.
            cluster_time (int): Number of clustering iterations.
            custom_label (str): Custom label for clustering.
            pca (str): PCA method to use.
            n_pcs (int): Number of principal components.
            pc_frac (float): Fraction of PCs to use.
            target_num (int): Number of clusters.
            cluster_thres (float): Threshold to filter small values in clustering.
            ot_cluster_thres (float): Threshold for OT clustering.
            emptygrid (np.array): Empty grid matrix.
            shape (str): Shape of spatial transcriptomics data ('square' or 'hexagon').
            device (str): Device for training ('cpu' or 'cuda').
        """
        self.scdata = scdata
        self.stdata = stdata
        self.bulkX = bulkX
        self.bulkY = bulkY
        self.genes = genes
        self.histology = histology
        self.cluster_time = cluster_time
        self.custom_label = custom_label
        self.pca = pca
        self.n_pcs = n_pcs
        self.pc_frac = pc_frac
        self.target_num = target_num
        self.cluster_thres = cluster_thres
        self.ot_cluster_thres = ot_cluster_thres
        self.emptygrid = emptygrid
        self.shape = shape
        self.device = device

        # Initialize matrices
        self.cluster_matrix = None
        self.spot_matrix = None

    def allocate(self,num_epochs=500):
        """
        Runs the Polyomino algorithm for mapping scRNA-seq data to spatial transcriptomics or bulk RNA datasets.
        Generates a matrix with probabilities of each cell belonging to specific clusters.
        """
        if self.stdata is not None:
            self.bulkX = generate_Xstrips(self.stdata)
            self.bulkY = generate_Ystrips(self.stdata)
            self.cluster_label = []

            if self.cluster_time > 0 and self.custom_label:
                self.cluster_label.append(self.custom_label)
                self.cluster_matrix = cluster_mapping(
                    self.scdata, self.stdata, self.genes, label=self.custom_label, device=self.device, thres=self.cluster_thres, num_epochs=num_epochs
                )
                torch.cuda.empty_cache()

                if self.cluster_thres is not None:
                    self.cluster_matrix[self.cluster_matrix < self.cluster_thres] = -1e15
            else:
                if self.cluster_time == 0:
                    print("Starting spot mapping.")
                    self.spot_matrix = spot_mapping(self.scdata, self.stdata, genes=self.genes, device=self.device, num_epochs=num_epochs)
                else:
                    self.cluster_matrix = 0

            if self.cluster_time > 0:
                for i in range(self.cluster_time):
                    label = f'clutimes_{i+1}'
                    self.cluster_label.append(label)
                    print(f"Running {i+1}-th clustering.")
                    random_cluster(
                        self.stdata, label, embed=self.pca, n_pcs=self.n_pcs, pc_frac=self.pc_frac,
                        samples_time=self.cluster_time, shape=self.shape, target_num=self.target_num
                    )
                    print(f"Running {i+1}-th mapping.")
                    self.cluster_matrix += cluster_mapping(
                        self.scdata, self.stdata, self.genes, label=label, device=self.device, thres=self.cluster_thres, num_epochs=num_epochs
                    )
                    torch.cuda.empty_cache()
                self.cluster_matrix /= self.cluster_time

                if self.cluster_thres is not None:
                    self.cluster_matrix[self.cluster_matrix < self.cluster_thres] = -1e15

                print("Starting spot mapping.")
                self.spot_matrix = spot_mapping(
                    self.scdata, self.stdata, self.cluster_matrix.values, genes=self.genes, device=self.device, num_epochs=num_epochs
                )
        else:
            print("Starting spot mapping.")
            self.spot_matrix = spot_mapping(
                self.scdata, None, self.cluster_matrix.values, self.bulkX, self.bulkY, genes=self.genes, device=self.device, num_epochs=num_epochs
            )

        self.scdata.obs['grid_prob'] = np.max(self.spot_matrix, axis=1)
        print("Spot mapping completed.")
        print("Mapping matrix saved in zm.spot_matrix")

    def ot_allocate(self):
        """
        Runs the optimal transport algorithm for mapping scRNA-seq data to spatial transcriptomics or bulk RNA datasets.
        """
        if self.stdata is not None:
            self.bulkX = sc.AnnData(X=generate_Xstrips(self.stdata), var=self.stdata.var)
            self.bulkY = sc.AnnData(X=generate_Ystrips(self.stdata), var=self.stdata.var)
            self.scLocX = strips_mapping(self.scdata, self.bulkX, self.genes, device=self.device, thres=self.ot_cluster_thres)
            self.scLocY = strips_mapping(self.scdata, self.bulkY, self.genes, device=self.device, thres=self.ot_cluster_thres)
            self.spot_matrix = solve_OT(self.emptygrid, self.scLocX, self.scLocY, thres=self.ot_cluster_thres, njob=8)
        else:
            self.scLocX = strips_mapping(self.scdata, self.bulkX, self.genes, device=self.device, thres=self.ot_cluster_thres)
            self.scLocY = strips_mapping(self.scdata, self.bulkY, self.genes, device=self.device, thres=self.ot_cluster_thres)
            self.spot_matrix = solve_OT(self.emptygrid, self.scLocX, self.scLocY, thres=self.ot_cluster_thres, njob=8)


class Polyomino3D:
    """
    This class runs the Polyomino algorithm on scRNA-seq data and Visium spatial transcriptomics datasets.
    It reconstructs new spatial transcriptomics data using the Polyomino algorithm.
    Attributes:
        scdata (anndata.AnnData): scRNA-seq data.
        stdata (anndata.AnnData): 3D spatial transcriptomics data.
        genes (list): List of genes to use for mapping.
        cluster_time (int): Number of clustering iterations.
        custom_label (str): Custom label for clustering.
        pca (str): PCA method to use.
        n_pcs (int): Number of principal components.
        pc_frac (float): Fraction of PCs to use.
        target_num (int): Number of clusters.
        cluster_thres (float): Threshold to filter small values in clustering.
        shape (str): Shape of the spatial transcriptomics data ('square' or 'hexagon').
        device (str): Device for training ('cpu' or 'cuda').
        cluster_matrix (np.ndarray): Matrix for clustering results.
        spot_matrix (np.ndarray): Matrix for spot mapping results.
    """

    def __init__(
        self, 
        scdata, 
        stdata,
        genes=None, 
        cluster_time=1, 
        custom_label=None,
        pca=None, 
        n_pcs=50, 
        pc_frac=0.5, 
        target_num=10,
        cluster_thres=None,
        shape="hexagon", 
        device='cpu'
    ):
        """
        Initializes the Polyomino3D class with provided parameters.

        Args:
            scdata (anndata.AnnData): scRNA-seq data.
            stdata (anndata.AnnData): 3D spatial transcriptomics data.
            genes (list): List of genes to use for mapping.
            cluster_time (int): Number of clustering iterations.
            custom_label (str): Custom label for clustering.
            pca (str): PCA method to use.
            n_pcs (int): Number of principal components.
            pc_frac (float): Fraction of PCs to use.
            target_num (int): Number of clusters.
            cluster_thres (float): Threshold to filter small values in clustering.
            shape (str): Shape of the spatial transcriptomics data ('square' or 'hexagon').
            device (str): Device for training ('cpu' or 'cuda').
        """
        self.scdata = scdata
        self.stdata = stdata
        self.genes = genes
        self.cluster_time = cluster_time
        self.custom_label = custom_label
        self.pca = pca
        self.n_pcs = n_pcs
        self.pc_frac = pc_frac
        self.target_num = target_num
        self.cluster_thres = cluster_thres
        self.shape = shape
        self.device = device

        # Initialize matrices
        self.cluster_matrix = None
        self.spot_matrix = None

    def allocate(self, num_epochs=500):
        """
        Runs the Polyomino algorithm for mapping scRNA-seq data to 3D spatial transcriptomics data.
        Generates a matrix with probabilities of each cell belonging to specific clusters.
        """
        self.bulkX = generate_Xstrips(self.stdata)
        self.bulkY = generate_Ystrips(self.stdata)
        self.bulkZ = generate_Zstrips(self.stdata)
        self.cluster_label = []

        if self.cluster_time > 0 and self.custom_label is not None:
            self.cluster_label.append(self.custom_label)
            self.cluster_matrix = cluster_mapping(
                self.scdata, self.stdata, self.genes, label=self.custom_label, device=self.device, thres=self.cluster_thres, num_epochs=num_epochs
            )
            torch.cuda.empty_cache()

            if self.cluster_thres is not None:
                self.cluster_matrix[self.cluster_matrix < self.cluster_thres] = -1e15
        else:
            if self.cluster_time == 0:
                print("Starting spot mapping.")
                self.spot_matrix = spot_mapping3D(self.scdata, self.stdata, genes=self.genes, device=self.device, num_epochs=num_epochs)
            else:
                self.cluster_matrix = 0

        if self.cluster_time > 0:
            for i in range(self.cluster_time):
                label = f'clutimes_{i}'
                self.cluster_label.append(label)
                print(f"Running {i}-th clustering.")
                random_cluster(
                    self.stdata, label, embed=self.pca, n_pcs=self.n_pcs, pc_frac=self.pc_frac,
                    samples_time=self.cluster_time, shape=self.shape, target_num=self.target_num
                )
                print(f"Running {i}-th mapping.")
                self.cluster_matrix += cluster_mapping(
                    self.scdata, self.stdata, self.genes, label=label, device=self.device, thres=self.cluster_thres, num_epochs=num_epochs
                )
                torch.cuda.empty_cache()

            self.cluster_matrix /= self.cluster_time

            if self.cluster_thres is not None:
                self.cluster_matrix[self.cluster_matrix < self.cluster_thres] = -1e15

            print("Starting spot mapping.")
            self.spot_matrix = spot_mapping3D(
                self.scdata, self.stdata, self.cluster_matrix.values, genes=self.genes, device=self.device, num_epochs=num_epochs
            )
        else:
            self.spot_matrix = spot_mapping3D(
                self.scdata, self.stdata, self.cluster_matrix.values, genes=self.genes, device=self.device, num_epochs=num_epochs
            )

        self.scdata.obs['grid_prob'] = np.max(self.spot_matrix, axis=1)
        print("Spot mapping completed.")
        print("Mapping matrix saved in zm.spot_matrix")


class Polyomino_custom:
    """
    This class runs a customized Polyomino algorithm on scRNA-seq data and spatial transcriptomics data.
    It uses a model-based approach to reconstruct single cell resolution spatial transcriptomics data.

    Attributes:
        scdata (anndata.AnnData): scRNA-seq data.
        stdata (anndata.AnnData): Spatial transcriptomics data.
        bulkX (anndata.AnnData): Bulk transcriptomics data (X strips).
        bulkY (anndata.AnnData): Bulk transcriptomics data (Y strips).
        genes (list): Genes to use for mapping.
        custom_label (str): Custom label for clustering.
        learning_rate (float): Learning rate for the optimization process.
        num_epochs (int): Number of epochs for training.
        lambda_gx1 (float): Regularization parameter for Gx1.
        lambda_gy1 (float): Regularization parameter for Gy1.
        lambda_r (float): Regularization parameter for the residual term (optional).
        device (str): Device for training ('cpu' or 'cuda').
        spot_matrix (np.ndarray): Matrix of spot mappings.
    """

    def __init__(
        self, 
        scdata, 
        stdata=None,
        genes=None, 
        custom_label=None,
        learning_rate=0.01,
        num_epochs=500,
        lambda_gx1=1,
        lambda_gy1=1,
        lambda_r=None,
        device='cpu'
    ):
        """
        Initializes the Polyomino_custom class with provided parameters.

        Args:
            scdata (anndata.AnnData): scRNA-seq data.
            stdata (anndata.AnnData): Spatial transcriptomics data.
            genes (list): List of genes to use for mapping.
            custom_label (str): Custom label for clustering.
            learning_rate (float): Learning rate for optimization.
            num_epochs (int): Number of epochs for training.
            lambda_gx1 (float): Regularization parameter for Gx1.
            lambda_gy1 (float): Regularization parameter for Gy1.
            lambda_r (float): Regularization parameter for the residual term (optional).
            device (str): Device for training ('cpu' or 'cuda').
        """
        self.scdata = scdata
        self.stdata = stdata
        self.bulkX = None
        self.bulkY = None
        self.genes = genes
        self.custom_label = custom_label
        self.spot_matrix = None
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.lambda_gx1 = lambda_gx1
        self.lambda_gy1 = lambda_gy1
        self.lambda_r = lambda_r
        self.device = device

    def allocate(self):
        """
        Runs the customized Polyomino algorithm on the scRNA-seq and spatial transcriptomics datasets.
        Produces a matrix with probabilities of each cell belonging to grids.
        """
        # Generate bulk data strips
        self.bulkX = sc.AnnData(X=generate_Xstrips(self.stdata), var=self.stdata.var)
        self.bulkY = sc.AnnData(X=generate_Ystrips(self.stdata), var=self.stdata.var)

        # Preprocess scRNA-seq data and get overlap genes
        preprocess(self.scdata, self.bulkX, self.genes)
        overlap_genes = self.scdata.uns["overlap"]

        # Extract data matrices
        S = np.array(self.scdata[:, overlap_genes].X.toarray(), dtype="float32")
        Gx1 = np.array(self.bulkX[:, overlap_genes].X, dtype="float32")
        Gx2 = np.array(self.bulkY[:, overlap_genes].X, dtype="float32")
        ST = self.stdata[:, overlap_genes]

        # Run enhanced cell-to-spot mapping
        self.spot_matrix = Cell2SpotsEnhance(
            S=S, ST=ST, Gx=Gx1, Gy=Gx2,
            lambda_gx1=self.lambda_gx1, lambda_gy1=self.lambda_gy1,
            custom_regions=self.custom_label, lambda_r=self.lambda_r,
            device=self.device
        ).fit(
            learning_rate=self.learning_rate, num_epochs=self.num_epochs, print_each=100
        )

        # Update scdata with the maximum spot probability
        self.scdata.obs['grid_prob'] = np.max(self.spot_matrix, axis=1)
        print("Finish spot mapping.")

def strips_mapping(
    scadata,
    bulk,
    genes=None,
    device='cpu',
    num_epochs=500,
    learning_rate=0.1,
    thres=None
):
    """
    Maps scRNA-seq data to bulk transcriptomics data using a model-based approach.

    Args:
        scadata (anndata.AnnData): scRNA-seq data.
        bulk (anndata.AnnData): Bulk transcriptomics data.
        genes (list, optional): List of genes to use for mapping. If None, all genes are used.
        device (str, optional): Device for training ('cpu' or 'cuda'). Default is 'cpu'.
        num_epochs (int, optional): Number of epochs for training. Default is 500.
        learning_rate (float, optional): Learning rate for training. Default is 0.1.
        thres (float, optional): Threshold to filter small values in the output matrix. If None, no filtering is applied.

    Returns:
        np.ndarray: Mapping matrix from scRNA-seq to bulk transcriptomics data.
    """
    # Preprocess data with the specified genes, if provided
    if genes is None:
        preprocess(scadata, bulk)
    else:
        preprocess(scadata, bulk, genes)
    
    # Get overlapping genes
    overlap_genes = scadata.uns["overlap"]
    
    # Convert scRNA-seq and bulk data to numpy arrays
    S = np.array(scadata[:, overlap_genes].X.toarray(), dtype="float32")
    G = np.array(bulk[:, overlap_genes].X, dtype="float32")
    
    # Initialize the cell-to-clusters mapper
    bulk_mapper = Cell2Clusters(scdata=S, clusters=G, device=device)
    
    # Fit the model to get the mapping matrix
    bulk_matrix = bulk_mapper.fit(num_epochs=num_epochs, learning_rate=learning_rate, print_each=100)
    
    # Apply threshold filtering if specified
    if thres is not None:
        bulk_matrix[bulk_matrix < thres] = 0
    return bulk_matrix

def cluster_mapping(
    scadata,
    stdata,
    genes=None,
    label='leiden',
    device='cpu',
    num_epochs=500,
    learning_rate=0.1,
    thres=None
):
    """
    Maps scRNA-seq data to specific regions based on spatial transcriptomics cluster labels.

    Args:
        scadata (anndata.AnnData): scRNA-seq data.
        stdata (anndata.AnnData): Spatial transcriptomics data.
        genes (list, optional): List of genes to use for mapping. If None, all genes are used.
        label (str): Column name in stdata.obs containing cluster labels.
        device (str, optional): Device for training ('cpu' or 'cuda'). Default is 'cpu'.
        num_epochs (int, optional): Number of epochs for training. Default is 500.
        learning_rate (float, optional): Learning rate for training. Default is 0.1.
        thres (float, optional): Threshold to filter small values in the output matrix. If None, no filtering is applied.

    Returns:
        np.ndarray: Mapping matrix from scRNA-seq to cluster transcriptomics data.
    """
    # Generate cluster expression data from spatial transcriptomics data
    cluster_bulk = generate_cluster_exp(stdata, label=label)
    cluster_bulk = sc.AnnData(X=cluster_bulk, var=stdata.var)
    
    # Preprocess the data with the specified genes, if provided
    if genes is None:
        preprocess(scadata, cluster_bulk)
    else:
        preprocess(scadata, cluster_bulk, genes)
    
    # Get overlapping genes
    overlap_genes = scadata.uns["overlap"]
    
    # Convert scRNA-seq and cluster data to numpy arrays
    scdata_arr = np.array(scadata[:, overlap_genes].X.toarray(), dtype="float32")
    clusters_arr = np.array(cluster_bulk[:, overlap_genes].X, dtype="float32")
    
    # Initialize the cell-to-clusters mapper
    cluster_mapper = Cell2Clusters(scdata=scdata_arr, clusters=clusters_arr, device=device)
    
    # Fit the model to get the mapping matrix
    cluster_mapper_matrix = cluster_mapper.fit(num_epochs=num_epochs, learning_rate=learning_rate, print_each=100)
    
    # Transform the cluster labels
    cell_clu = one_hot_encode_labels(stdata.obs[label])
    
    # Compute the mapping matrix
    cell_spot = cluster_mapper_matrix @ cell_clu.T
    
    # Apply threshold filtering if specified
    if thres is not None:
        return (cell_spot > thres).astype(int)
    
    return cell_spot

def spot_mapping(
    scadata,
    stdata,
    cluster_mapping_matrix=None,
    x_bulk=None,
    y_bulk=None,
    genes=None,
    device='cpu',
    num_epochs=500,
    learning_rate=0.1
):
    """
    Maps single-cell RNA-seq data to spatial spots using spatial transcriptomics data.

    Args:
        scadata (anndata.AnnData): scRNA-seq data.
        stdata (anndata.AnnData): Spatial transcriptomics data.
        cluster_mapping_matrix (np.ndarray, optional): Matrix mapping scRNA-seq data to clusters.
        x_bulk (anndata.AnnData, optional): Bulk transcriptomics data in X direction. If None, generated from stdata.
        y_bulk (anndata.AnnData, optional): Bulk transcriptomics data in Y direction. If None, generated from stdata.
        genes (list, optional): List of genes to use for mapping. If None, all overlapping genes are used.
        device (str, optional): Device for training ('cpu' or 'cuda'). Default is 'cpu'.
        num_epochs (int, optional): Number of epochs for training. Default is 500.
        learning_rate (float, optional): Learning rate for training. Default is 0.1.

    Returns:
        np.ndarray: Mapping matrix from single cells to spatial spots.
    """
    # Generate bulk data if not provided
    if x_bulk is None and y_bulk is None:
        x_bulk = sc.AnnData(X=generate_Xstrips(stdata), var=stdata.var)
        y_bulk = sc.AnnData(X=generate_Ystrips(stdata), var=stdata.var)
    
    # Preprocess the scRNA-seq data
    preprocess(scadata, x_bulk, genes)
    
    # Get overlapping genes
    overlap_genes = scadata.uns["overlap"]
    
    # Convert data to numpy arrays
    S = np.array(scadata[:, overlap_genes].X.toarray(), dtype="float32")
    Gx1 = np.array(x_bulk[:, overlap_genes].X, dtype="float32")
    Gx2 = np.array(y_bulk[:, overlap_genes].X, dtype="float32")
    
    # Perform mapping to spatial spots
    if stdata is None:
        # Use cell2spots_strips if spatial data is not provided
        mapping_matrix = Cell2SpotsStrips(
            S=S,
            Gx=Gx1,
            Gy=Gx2,
            cluster_matrix=cluster_mapping_matrix,
            device=device
        ).fit(learning_rate=learning_rate, num_epochs=num_epochs, print_each=100)
    else:
        # Use cell2spots if spatial data is provided
        ST = stdata[:, overlap_genes]
        mapping_matrix = Cell2Spots(
            S=S,
            ST=ST,
            Gx=Gx1,
            Gy=Gx2,
            cluster_matrix=cluster_mapping_matrix,
            device=device
        ).fit(learning_rate=learning_rate, num_epochs=num_epochs, print_each=100)
    
    return mapping_matrix

def spot_mapping3D(
    scadata,
    stdata,
    cluster_mapping_matrix=None,
    genes=None,
    device='cpu',
    num_epochs=500,
    learning_rate=0.1
):
    """
    Maps single-cell RNA-seq data to spatial spots in a 3D spatial transcriptomics setup.

    Args:
        scadata (anndata.AnnData): scRNA-seq data.
        stdata (anndata.AnnData): Spatial transcriptomics data.
        cluster_mapping_matrix (np.ndarray, optional): Matrix mapping scRNA-seq data to clusters.
        genes (list, optional): List of genes to use for mapping. If None, all overlapping genes are used.
        device (str, optional): Device for training ('cpu' or 'cuda'). Default is 'cpu'.
        num_epochs (int, optional): Number of epochs for training. Default is 500.
        learning_rate (float, optional): Learning rate for training. Default is 0.1.

    Returns:
        np.ndarray: Mapping matrix from single cells to spatial spots.
    """
    # Generate bulk data for X, Y, and Z directions
    x_bulk = generate_Xstrips(stdata)
    y_bulk = generate_Ystrips(stdata)
    z_bulk = generate_Zstrips(stdata)
    
    # Create AnnData objects for bulk data
    bx1 = sc.AnnData(X=x_bulk, var=stdata.var)
    bx2 = sc.AnnData(X=y_bulk, var=stdata.var)
    bx3 = sc.AnnData(X=z_bulk, var=stdata.var)
    
    # Preprocess scRNA-seq data
    preprocess(scadata, bx1, genes)
    
    # Get overlapping genes
    overlap_genes = scadata.uns["overlap"]
    
    # Convert data to numpy arrays
    S = np.array(scadata[:, overlap_genes].X.toarray(), dtype="float32")
    Gx1 = np.array(bx1[:, overlap_genes].X, dtype="float32")
    Gx2 = np.array(bx2[:, overlap_genes].X, dtype="float32")
    Gx3 = np.array(bx3[:, overlap_genes].X, dtype="float32")
    
    # Extract spatial transcriptomics data
    ST = stdata[:, overlap_genes]
    
    # Perform mapping to spatial spots
    mapping_matrix = Cell2VisiumSpots(
        S=S,
        ST=ST,
        Gx=Gx1,
        Gy=Gx2,
        Gz=Gx3,
        cluster_matrix=cluster_mapping_matrix,
        device=device
    ).fit(learning_rate=learning_rate, num_epochs=num_epochs, print_each=100)
    return mapping_matrix

def sc2sc(
    scadata,
    stdata_raw,
    mapping_matrix,
    thres=0.1,
    method='max'
):
    """
    Maps single-cell RNA-seq data to spatial cells within defined grids.

    Args:
        scadata (anndata.AnnData): scRNA-seq data.
        stdata_raw (anndata.AnnData): Spatial transcriptomics data with single cell resolution.
        mapping_matrix (np.ndarray): Matrix mapping single cells to spatial spots.
        thres (float, optional): Threshold to filter small values. Default is 0.1.
        method (str, optional): Method to select the cells. Options are 'max' and 'lap'. Default is 'max'.

    Returns:
        anndata.AnnData: scRNA-seq data with spatial information, including similarity scores.
    """
    st_x = []
    st_y = []
    select_ct_index = []
    select_gep = []
    similarity_scores = [] 

    # Determine number of cells to select for each grid
    select_cells_num = np.sum(mapping_matrix > thres, axis=0)
    raw_spot_index = stdata_raw.obs['spot_index'].unique()
    overlap_genes = scadata.uns["overlap"]

    for i in tqdm.tqdm(raw_spot_index):
        # Find index of the current spot
        sorted_spot_indices = np.argsort(raw_spot_index)
        index_of_spot = np.where(raw_spot_index[sorted_spot_indices] == i)[0][0]

        # Extract relevant data for the current spot
        st_temp = stdata_raw[stdata_raw.obs['spot_index'] == i][:, overlap_genes]
        num_cells = select_cells_num[index_of_spot]
        if num_cells == 0:
            num_cells = st_temp.shape[0]

        # Select the top cells based on mapping matrix
        sc_temp = scadata[np.argsort(mapping_matrix[:, index_of_spot], axis=0)[-num_cells:], :][:, overlap_genes]
        cos_result = cos_s(sc_temp.X, st_temp.X).T

        # Select cells based on the chosen method
        if method == 'max':
            select_index = np.argmax(cos_result, axis=1)
            similarity_scores.extend(np.max(cos_result, axis=1).tolist())
        elif method == 'lap':
            import lap
            select_index = lap.lapjv(1 - cos_result, extend_cost=True)[1]
            similarity_scores.extend([cos_result[row, col] for row, col in enumerate(select_index)])
        else:
            raise ValueError("Method must be 'max' or 'lap'.")

        # Collect data for the selected cells
        sc_select = scadata[np.argsort(mapping_matrix[:,index_of_spot], axis=0)[-num_cells:],:][select_index]
        select_ct_index.extend(sc_select.obs.index.tolist())
        select_gep.append(sc_select.X.toarray())
        st_x.extend(st_temp.obs.x.values.tolist())
        st_y.extend(st_temp.obs.y.values.tolist())

    # Create an AnnData object for the result
    cell_alocated_data = sc.AnnData(
        np.vstack(select_gep),
        obs=pd.DataFrame(index=select_ct_index),
        var=scadata.var
    )
    cell_alocated_data.obsm['spatial'] = np.array([st_x, st_y]).T
    cell_alocated_data.obs = scadata[cell_alocated_data.obs.index, :].obs
    cell_alocated_data.obsm['similarity_score'] = np.array(similarity_scores)

    return cell_alocated_data

def sc2st(
    scadata,
    stdata,
    stdata_raw,
    mapping_matrix,
    thres=0.01):
    """
    Mapping the single cells to the single cells of maximized proability spatial spots.
    Args:
        scadata (anndata.AnnData): scRNA-seq data
        stdata_raw (anndata.AnnData): spatial transcriptomics data with single cell resolution.
        mapping_matrix (np.array): sc to spots mapping matrix
        thres (floot): thres to filter small values
        method (str): method to select the cells
    Return:
        cell_alocated_data (anndata.AnnData): scRNA-seq data with spatial information
    """
    overlap_genes = scadata.uns["overlap"]
    scadata = scadata[:,overlap_genes].copy()
    stdata = stdata[:,overlap_genes].copy()
    stdata_raw = stdata_raw[:,overlap_genes].copy()
    stdata_raw.obs['spot_index'] = stdata_raw.obs['spot_index'].astype('str')
    spot_index = stdata[np.argmax(mapping_matrix, axis=1)].obs.index
    raw_spot_index = stdata_raw.obs['spot_index'].unique()
    sc_location = []
    sc_index = []
    import lap
    for i in tqdm.tqdm(raw_spot_index):
        st_temp = stdata_raw[stdata_raw.obs['spot_index'] == i]
        sc_temp = scadata[spot_index == i]
        if len(sc_temp) == 0 :
            continue
        cos_result = cos_s(sc_temp.X,st_temp.X)
        select_index = lap.lapjv(1-cos_result,extend_cost=True)[1]
        sc_index.extend(sc_temp.obs.index.tolist())
        sc_location.extend(st_temp[select_index].obsm['spatial'])
    cell_alocated_data = scadata[sc_index]
    cell_alocated_data.obsm['spatial'] = np.array(sc_location)
    return cell_alocated_data

def main_Polyomino():
    args = argument_parser()
    sc_adata = sc.read_h5ad(args.sc)
    st_adata = sc.read_h5ad(args.st)
    stdata_spot = generate_grid(st_adata, width=args.w)
    if args.custom_region is not None:
        zm_obj = Polyomino(sc_adata, stdata_spot, 
                     cluster_time=args.cluster_time, cluster_thres=args.cluster_thres,custom_label = args.custom_region,
                     device=args.device)
    else:
        zm_obj = Polyomino(sc_adata, stdata_spot, 
                     cluster_time=args.cluster_time, cluster_thres=args.cluster_thres,
                     device=args.device)
    zm_obj.allocate()
    cell_allocated_data = sc2sc(sc_adata, st_adata, zm_obj.spot_matrix, 
                                        thres=args.thres if args.thres else 0.1, 
                                        method=args.method)
    cell_allocated_data.write_h5ad(args.output)
if __name__ == "__main__":
    main_Polyomino()
