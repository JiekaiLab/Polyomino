import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from sklearn.decomposition import PCA
from scipy import sparse
from scipy.sparse import issparse
import random
from math import sin, cos, pi, sqrt
import matplotlib.pyplot as plt
import seaborn as sns
from ..extension.custom_SpaGCN import *

def rotate_coordinates(Mx, My, angle):
    """
    Rotate coordinates by a given angle.

    Args:
        Mx (np.ndarray): X-coordinates.
        My (np.ndarray): Y-coordinates.
        angle (float): Angle of rotation in radians.

    Returns:
        tuple: Rotated X and Y coordinates.
    """
    Mx_rotated = Mx * cos(angle) - My * sin(angle)
    My_rotated = Mx * sin(angle) + My * cos(angle)
    return Mx_rotated, My_rotated

def one_hot_encode_labels(labels: pd.Series):
    """
    Convert categorical labels to one-hot encoded format.

    Args:
        labels (pd.Series): Series containing categorical labels.

    Returns:
        pd.DataFrame: One-hot encoded labels.
    """
    return pd.get_dummies(labels)

def transfer_labels(mapping_matrix: np.ndarray, stdata: anndata.AnnData, labels: pd.Series):
    """
    Transfer labels from single-cell data to spatial transcriptomics data.

    Args:
        mapping_matrix (np.ndarray): Matrix mapping cells to spots.
        stdata (anndata.AnnData): Spatial transcriptomics data.
        labels (pd.Series): Cell labels.

    Returns:
        pd.DataFrame: Probability distribution of cell types for each spot.
    """
    one_hot_labels = one_hot_encode_labels(labels)
    cell_type_prob = pd.DataFrame(mapping_matrix.T @ one_hot_labels)
    cell_type_prob.index = stdata.obs.index
    return cell_type_prob

def create_spot_anndata(scdata, stdata, mapping_matrix, sclabel):
    """
    Create an AnnData object at spot resolution.

    Args:
        scdata (anndata.AnnData): Single-cell data.
        stdata (anndata.AnnData): Spatial transcriptomics data.
        mapping_matrix (np.ndarray): Matrix mapping cells to spots.
        sclabel (str): Column name of cell labels in scdata.

    Returns:
        anndata.AnnData: Spot resolution AnnData object.
    """
    spot_data = sc.AnnData(X=mapping_matrix.T @ scdata.X, var=scdata.var, obsm=stdata.obsm)
    spot_data.obsm['celltype_matrix'] = transfer_labels(mapping_matrix, stdata, scdata.obs[sclabel])
    spot_data.obsm['mapping_matrix'] = sparse.csr_matrix(mapping_matrix.T, dtype=float)
    return spot_data

def calculate_optimal_width(stdata_raw, target_spots=10000):
    """
    Calculate the optimal grid width to achieve a target number of spots.

    Args:
        stdata_raw (anndata.AnnData): Raw spatial data with single-cell resolution.
        target_spots (int): Desired number of spots.

    Returns:
        int: Optimal grid width.
    """
    spatial_coords = stdata_raw.obsm['spatial']
    stdata_raw_x = spatial_coords[:, 0] - spatial_coords[:, 0].min()
    stdata_raw_y = spatial_coords[:, 1] - spatial_coords[:, 1].min()

    x_range = stdata_raw_x.max() - stdata_raw_x.min()
    y_range = stdata_raw_y.max() - stdata_raw_y.min()

    total_area = x_range * y_range
    optimal_width = np.sqrt(total_area / target_spots)

    return int(optimal_width)

def generate_grid(stdata_raw, width=None):
    """
    Generate spatial transcriptomics data at specific grid resolution.

    Args:
        stdata_raw (anndata.AnnData): Raw spatial data with single-cell resolution.
        width (int): Width of the spots.

    Returns:
        anndata.AnnData: Spot resolution spatial transcriptomics data.
    """
    spatial_coords = stdata_raw.obsm['spatial']
    stdata_raw.obs['x'] = spatial_coords[:, 0]
    stdata_raw.obs['y'] = spatial_coords[:, 1]
    
    if width is None:
        width = calculate_optimal_width(stdata_raw)

    stdata_raw_x = stdata_raw.obs['x'] - stdata_raw.obs['x'].min()
    stdata_raw_y = stdata_raw.obs['y'] - stdata_raw.obs['y'].min()

    nx = int((stdata_raw_x.max() - stdata_raw_x.min()) // width) + 1
    ny = int((stdata_raw_y.max() - stdata_raw_y.min()) // width) + 1

    stdata_raw.obs['spot_index'] = (np.floor(stdata_raw_x / width) + np.floor(stdata_raw_y / width) * nx).astype(int)
    stdata_raw_df = stdata_raw.to_df()
    stdata_raw_df['spot_index'] = stdata_raw.obs['spot_index'].values
    stdata_raw_ep = stdata_raw_df.groupby('spot_index').sum()

    stdata = sc.AnnData(X=stdata_raw_ep.values, var=stdata_raw.var)
    stdata.obs.index = stdata_raw_ep.index.astype(str)
    stdata.obs['x'] = stdata_raw_ep.index % nx
    stdata.obs['y'] = stdata_raw_ep.index // nx
    stdata.obsm['spatial'] = stdata.obs[['x', 'y']].values

    stdata.obs['array_col'] = stdata.obsm['spatial'][:, 0].astype('int32')
    stdata.obs['array_row'] = stdata.obsm['spatial'][:, 1].astype('int32')
    stat_count = pd.value_counts(stdata_raw.obs['spot_index'])
    median_num = np.median(pd.value_counts(stdata_raw.obs['spot_index']).values)
    sns.histplot(pd.value_counts(stdata_raw.obs['spot_index']))
    plt.xlabel('Cell numbers of grid')
    plt.title('Median of cells in one grid: ' + str(median_num))
    plt.grid(None)
    return stdata

def generate_spot_matrix(stdata_raw, width):
    """
    Generate the spatial spot matrix representing cell presence in spots.

    Args:
        stdata_raw (anndata.AnnData): Raw spatial data with single-cell resolution.
        width (int): Width of the spots.

    Returns:
        np.ndarray: 3D matrix representing spatial spots.
    """
    spatial_coords = stdata_raw.obsm['spatial']
    stdata_raw.obs['x'] = spatial_coords[:, 0] - spatial_coords[:, 0].min()
    stdata_raw.obs['y'] = spatial_coords[:, 1] - spatial_coords[:, 1].min()

    nx = int((stdata_raw.obs['x'].max() - stdata_raw.obs['x'].min()) // width) + 1
    ny = int((stdata_raw.obs['y'].max() - stdata_raw.obs['y'].min()) // width) + 1

    spot_matrix = np.zeros((len(stdata_raw), nx * ny), dtype=int)

    for i, cell in enumerate(stdata_raw):
        x, y = cell.obs['x'].values, cell.obs['y'].values
        spot_idx = int(np.floor(x / width) + np.floor(y / width) * nx)
        spot_matrix[i, spot_idx] = 1

    spot_matrix_reshaped = spot_matrix.reshape((len(stdata_raw), ny, nx))
    return spot_matrix_reshaped[:, ::-1, :]  # Reverse the y-axis to match the typical plot orientation

def generate_empty_grid(stdata_raw, stdata):
    """
    Generate an empty grid matrix based on spatial data.

    Args:
        stdata_raw (anndata.AnnData): Raw spatial data with single-cell resolution.
        stdata (anndata.AnnData): Spot resolution spatial data.

    Returns:
        np.ndarray: Empty grid matrix.
    """
    grid_shape = (int(stdata.obs['y'].max()) + 1, int(stdata.obs['x'].max()) + 1)
    empty_grid = np.zeros(grid_shape, dtype=int)

    for idx in stdata_raw.obs['spot_index'].unique():
        empty_grid.flat[int(idx)] = 1

    return empty_grid

def generate_Xstrips(stdata):
    """this function is to generate X strips, which is used for calculating the correlation matrix
    Args:
        stdata (anndata.AnnData): Visium data
    Returns:
        the function will return the X strips and update the 'x' column in stdata.obs
    """
    slides_x =  []
    x_min = stdata.obs['array_col'].min()
    x_max = stdata.obs['array_col'].max()
    lenth = x_max-x_min+1
    j = 0
    for i in range(np.int32(lenth)):
        idx = abs(stdata.obs['array_col'] - x_min - i) < 0.1
        if sum(idx)==0:
            continue
        slides_x.append(np.asarray(stdata.X[idx,:].sum(axis =0)).reshape(-1)) 
        stdata.obs.loc[idx,'x'] = j
        j+=1
    return np.array(slides_x)

def generate_Ystrips(stdata):
    """this function is to generate Y strips, which is used for calculating the correlation matrix
    Args:
        stdata (anndata.AnnData): Visium data
    Returns:
        the function will return the Y strips and update the 'y' column in stdata.obs
    """
    slides_y =  []
    y_min = stdata.obs['array_row'].min()
    y_max = stdata.obs['array_row'].max()
    lenth = y_max-y_min+1
    j = 0
    for i in range(np.int32(lenth)):
        idx = abs(stdata.obs['array_row'] - y_min - i) < 0.1
        if sum(idx)==0:
            continue
        slides_y.append(np.asarray(stdata.X[idx,:].sum(axis =0)).reshape(-1)) 
        stdata.obs.loc[idx,'y'] = j
        j+=1
    return np.array(slides_y)

def generate_Zstrips(stdata):
    """Generate X strips
    Args:
        stdata (anndata): Visium data
    Returns:
        BYn: _description_
    """
    slides_z =  []
    # Mx_p,My_p = rotation_exp(stdata.obs['array_col'],stdata.obs['array_row'],angle=pi/4)
    Mx_p,My_p = rotation_exp(stdata.obs['array_col'],stdata.obs['array_row'],angle=pi/4)
    lenth = np.round((Mx_p.max()-Mx_p.min())/sqrt(2))+1
    start_value = Mx_p.min()
    j = 0 
    for i in range(np.int32(lenth)):
        idx = abs(Mx_p - start_value - i * sqrt(2)) < 0.1
        if sum(idx)==0:
            continue
        slides_z.append(np.asarray(stdata.X[idx,:].sum(axis =0)).reshape(-1)) 
        stdata.obs.loc[idx,'z'] = j
        j += 1
    return np.array(slides_z)

def generate_cluster_exp(stdata, label):
    """
    Generate a cluster expression matrix.

    Args:
        stdata (anndata.AnnData): Spatial matrix containing the data.
        label (str): The column name of the clusters label stored in the spatial data's `.obs`.

    Returns:
        np.ndarray: A matrix where each row corresponds to the sum of expression values for each cluster.
    """
    clusters = []
    for cluster_label in stdata.obs[label].unique():
        cluster_data = stdata[stdata.obs[label] == cluster_label]
        clusters.append(np.asarray(cluster_data.X.sum(axis=0)).reshape(-1))
    return np.array(clusters)


def preprocess(scdata, stdata, genes=None):
    """
    Preprocess the single-cell and spatial data by filtering out genes not expressed in any cells and identifying overlapping genes.

    Args:
        scdata (anndata.AnnData): Single-cell RNA-seq data.
        stdata (anndata.AnnData): Spatial transcriptomics data.
        genes (list, optional): List of genes to consider. If None, all genes are considered.

    Returns:
        None: The function updates `scdata.uns["overlap_genes"]` with the overlapping genes.
    """
    sc.pp.filter_genes(scdata, min_cells=1)
    if genes is None:
        genes = scdata.var.index
    overlap_genes = list(set(genes) & set(scdata.var.index) & set(stdata.var.index))
    scdata.uns["overlap"] = overlap_genes


def random_cluster(stdata, label, embed=None, n_pcs=50, pc_frac=0.5, samples_time=1, shape="square", target_num=10):
    """
    Generate random cluster labels using PCA and SpaGCN clustering.

    Args:
        stdata (anndata.AnnData): Spatial matrix containing the data.
        label (str): Column name to store the generated cluster labels in `stdata.obs`.
        embed (np.ndarray, optional): Precomputed embeddings. If None, PCA will be performed.
        n_pcs (int, optional): Number of principal components to use for clustering. Default is 50.
        pc_frac (float, optional): Fraction of PCs to randomly select. Default is 0.5.
        samples_time (int, optional): Number of times to sample. Default is 1.
        shape (str, optional): Shape of the spot. Default is "square". Other options are "hex" for Visium data.
        target_num (int, optional): Target number of clusters. Default is 10.

    Returns:
        None: The function updates `stdata.obs` with the generated cluster labels.
    """
    if embed is None:
        pca = PCA(n_components=n_pcs)
        if issparse(stdata.X):
            pca.fit(stdata.X.A)
            embed = pca.transform(stdata.X.A)
        else:
            pca.fit(stdata.X)
            embed = pca.transform(stdata.X)
    
    selected_pcs_indices = sorted(random.sample(range(n_pcs), int(n_pcs * pc_frac)))
    selected_pcs = embed[:, selected_pcs_indices]
    SpaGCN_cluster(selected_pcs, stdata, label, shape=shape, target_num=target_num)


def SpaGCN_cluster(selected_pcs, stdata, label, shape="square", target_num=10):
    """
    Perform clustering using SpaGCN on the selected principal components.

    Args:
        selected_pcs (np.ndarray): Selected principal components for clustering.
        stdata (anndata.AnnData): Spatial matrix containing the data.
        label (str): Column name to store the generated cluster labels in `stdata.obs`.
        shape (str, optional): Shape of the spot. Default is "square". Other options are "hex" for Visium data.
        target_num (int, optional): Target number of clusters. Default is 10.

    Returns:
        None: The function updates `stdata.obs` with the refined cluster labels.
    """
    prefilter_specialgenes(stdata)

    x_pixel = stdata.obs["x"].tolist()
    y_pixel = stdata.obs["y"].tolist()
    adj_no_img = calculate_adj_matrix(x=x_pixel, y=y_pixel, histology=False)
    
    l = search_l(p=0.5, adj=adj_no_img, start=0.01, end=1000, tol=0.01, max_run=100)
    
    r_seed = t_seed = n_seed = 100
    res = search_res(
        stdata, adj_no_img, selected_pcs, l,
        target_num=target_num, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20,
        r_seed=r_seed, t_seed=t_seed, n_seed=n_seed
    )
    
    clf = SpaGCN()
    clf.set_l(l)
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    
    clf.train(
        stdata, adj_no_img, selected_pcs,
        init_spa=True, init="leiden", res=res, tol=5e-3, lr=0.05, max_epochs=200
    )
    y_pred, prob = clf.predict()
    stdata.obs[label]= y_pred
    stdata.obs[label]=stdata.obs[label].astype('category')
    refined_pred=refine(sample_id=stdata.obs.index.tolist(), pred=stdata.obs[label].tolist(), dis=adj_no_img, shape=shape)
    stdata.obs[label]=refined_pred
    stdata.obs[label]=stdata.obs[label].astype('category')
    
def sc_accu(recon_data, raw_st, recon_label, raw_label, nn=5, thres=1):
    """
    Calculate the accuracy of the reconstructed cell types by comparing them to the raw spatial data.

    Args:
        recon_data (anndata.AnnData): Reconstructed spatial data.
        raw_st (anndata.AnnData): Raw spatial data.
        recon_label (str): Column name for the reconstructed cell type labels in `recon_data.obs`.
        raw_label (str): Column name for the true cell type labels in `raw_st.obs`.
        nn (int, optional): Number of nearest neighbors to consider. Default is 5.
        thres (int, optional): Threshold for the number of matches required to consider a cell type correct. Default is 1.

    Returns:
        pd.DataFrame: A DataFrame containing the neighbor analysis and match counts.
    """
    from sklearn.neighbors import KDTree

    # Create KDTree for nearest neighbor search
    kdt = KDTree(recon_data.obsm['spatial'], leaf_size=30, metric='euclidean')
    nn_result = kdt.query(recon_data.obsm['spatial'], k=nn, return_distance=False)

    # Create a mapping from the raw spatial data indices to cell type labels
    mapping_dict = dict(zip(range(len(raw_st)), raw_st.obs[raw_label]))

    # Create a DataFrame with neighbor indices and map them to true cell types
    column_names = [f"neighbor_{i}" for i in range(nn)]
    results = pd.DataFrame(nn_result, columns=column_names).applymap(mapping_dict.get)

    # Add reconstructed and true cell type labels to the DataFrame
    results['Recon_ct'] = recon_data.obs[recon_label].values
    results['True_ct'] = raw_st.obs[raw_label].values

    # Calculate the number of matches with the reconstructed cell type
    results['Matches'] = results[column_names].eq(results['Recon_ct'], axis=0).sum(axis=1)

    # Determine if the number of matches meets or exceeds the threshold
    results['Values_Equal'] = results['Matches'] >= thres

    return results


def fast_pearson(v1, v2):
    """
    Compute Pearson correlation coefficients quickly between two sets of vectors.

    Args:
        v1 (np.ndarray): First set of vectors, shape (n_samples, n_features).
        v2 (np.ndarray): Second set of vectors, shape (n_samples, n_features).

    Returns:
        np.ndarray: Pearson correlation coefficients matrix.
    """
    n = v1.shape[0]
    sums = np.multiply.outer(v2.sum(axis=0), v1.sum(axis=0))
    stds = np.multiply.outer(v2.std(axis=0), v1.std(axis=0))
    correlation = (v2.T.dot(v1) - sums / n) / stds / n

    return correlation


def fast_spearman(v1, v2):
    """
    Compute Spearman correlation coefficients quickly between two sets of vectors.

    Args:
        v1 (np.ndarray): First set of vectors, shape (n_samples, n_features).
        v2 (np.ndarray): Second set of vectors, shape (n_samples, n_features).

    Returns:
        np.ndarray: Spearman correlation coefficients matrix.
    """
    v1 = pd.DataFrame(v1).rank().values
    v2 = pd.DataFrame(v2).rank().values
    n = v1.shape[0]
    sums = np.multiply.outer(v2.sum(axis=0), v1.sum(axis=0))
    stds = np.multiply.outer(v2.std(axis=0), v1.std(axis=0))
    correlation = (v2.T.dot(v1) - sums / n) / stds / n

    return correlation
