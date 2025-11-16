import numpy as np
import pandas as pd
import torch
from torch.nn.functional import softmax, cosine_similarity

class Cell2Clusters:
    """
    Maps each cell to its region based on pseudo spatial expression matrix and reconstructed cell expression vectors.
    
    Attributes:
        scdata (torch.Tensor): Single-cell expression matrix.
        clusters (torch.Tensor): Spatial cluster expression matrix.
        lambda_g1 (float): Weight for the first term in the loss function.
        lambda_g2 (float): Weight for the second term in the loss function.
        device (str): Device to use ('cpu' or 'cuda').
        random_state (int): Random seed for initialization.
    """
    
    def __init__(
        self,
        scdata,
        clusters,
        lambda_g1=1.0,
        lambda_g2=1.0,
        device="cpu",
        random_state=1597,
    ):
        """
        Initializes the Cell2Clusters instance.
        
        Args:
            scdata (ndarray): Single-cell expression matrix.
            clusters (ndarray): Regional spatial cluster expression matrix.
            lambda_g1 (float): Weight for the first term in the loss function.
            lambda_g2 (float): Weight for the second term in the loss function.
            device (str): Device to use ('cpu' or 'cuda').
            random_state (int): Random seed for initialization.
        """
        self.scdata = torch.tensor(scdata, device=device, dtype=torch.float32)
        self.clusters = torch.tensor(clusters, device=device, dtype=torch.float32)
        self.lambda_g1 = lambda_g1
        self.lambda_g2 = lambda_g2
        self.device = device
        self.random_state = random_state
        np.random.seed(seed=self.random_state)
        self.cluster_matrix = torch.tensor(
            np.random.normal(0, 1, (scdata.shape[0], clusters.shape[0])),
            device=device, requires_grad=True, dtype=torch.float32
        )

    def _loss_fn(self, verbose=False):
        """
        Computes the loss function for the cell-to-cluster mapping.
        
        Args:
            verbose (bool): If True, prints the loss.
        
        Returns:
            tuple: Total loss, cluster cosine similarity (term 1), and cluster cosine similarity (term 2).
        """
        cluster_probs = softmax(self.cluster_matrix, dim=1)
        G_rec = torch.matmul(cluster_probs.t(), self.scdata)
        cluster_cos1 = self.lambda_g1 * cosine_similarity(G_rec, self.clusters, dim=0).mean()
        cluster_cos2 = self.lambda_g2 * cosine_similarity(G_rec, self.clusters, dim=1).mean()
        expression_term = cluster_cos1 + cluster_cos2
        total_loss = -expression_term
        
        if verbose:
            print(f"Total loss: {total_loss.item():.3f}")
        
        return total_loss, cluster_cos1, cluster_cos2

    def fit(self, num_epochs = 500,learning_rate=0.1, print_each=100):
        """
        Runs the optimization process and returns the mapping matrix.
        
        Args:
            num_epochs (int): Number of epochs to train.
            learning_rate (float): Learning rate for the optimizer.
            print_each (int): Frequency of loss printing.
        
        Returns:
            ndarray: Optimized cluster mapping matrix (number_cells x number_clusters).
        """
        torch.manual_seed(seed=self.random_state)
        optimizer = torch.optim.Adam([self.cluster_matrix], lr=learning_rate)
        
        for t in range(num_epochs):
            loss, _, _ = self._loss_fn(verbose=(print_each is not None and t % print_each == 0))
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        
        with torch.no_grad():
            cluster_probs = softmax(self.cluster_matrix, dim=1).cpu().numpy()
        
        # Clear CUDA memory if using GPU
        if self.device == 'cuda':
            torch.cuda.empty_cache()
        
        return cluster_probs

class Cell2Spots:
    """
    Maps each cell to its spot based on multi layer regional spatial constraints and reconstructed cell expression vectors.
    
    Attributes:
        S (torch.Tensor): Single-cell expression matrix.
        ST (AnnData): Spatial data.
        Gx (torch.Tensor): X-axis spatial expression matrix.
        Gy (torch.Tensor): Y-axis spatial expression matrix.
        lambda_gx1 (float): Weight for the first term in X-axis loss function.
        lambda_gx2 (float): Weight for the second term in X-axis loss function.
        lambda_gy1 (float): Weight for the first term in Y-axis loss function.
        lambda_gy2 (float): Weight for the second term in Y-axis loss function.
        cluster_matrix (torch.Tensor): Initial cluster matrix for optimization.
        device (str): Device to use ('cpu' or 'cuda').
        random_state (int): Random seed for initialization.
    """
    
    def __init__(
        self,
        S,
        ST,
        Gx,
        Gy,
        lambda_gx1=1,
        lambda_gx2=0,
        lambda_gy1=1,
        lambda_gy2=0,
        cluster_matrix=None,
        device="cpu",
        random_state=1597,
    ):
        """
        Initializes the Cell2Spots instance.
        
        Args:
            S (ndarray): Single-cell expression matrix.
            ST (AnnData): Spatial data.
            Gx (ndarray): X-axis spatial expression matrix.
            Gy (ndarray): Y-axis spatial expression matrix.
            lambda_gx1 (float): Weight for the first term in X-axis loss function.
            lambda_gx2 (float): Weight for the second term in X-axis loss function.
            lambda_gy1 (float): Weight for the first term in Y-axis loss function.
            lambda_gy2 (float): Weight for the second term in Y-axis loss function.
            cluster_matrix (ndarray or None): Initial cluster matrix for optimization.
            device (str): Device to use ('cpu' or 'cuda').
            random_state (int): Random seed for initialization.
        """
        self.S = torch.tensor(S, device=device, dtype=torch.float32)
        self.ST = ST
        self.Gx = torch.tensor(Gx, device=device, dtype=torch.float32)
        self.Gy = torch.tensor(Gy, device=device, dtype=torch.float32)
        self.lambda_gx1 = lambda_gx1
        self.lambda_gx2 = lambda_gx2
        self.lambda_gy1 = lambda_gy1
        self.lambda_gy2 = lambda_gy2
        self.random_state = random_state
        np.random.seed(seed=self.random_state)
        
        if cluster_matrix is None:
            self.cluster_matrix = torch.tensor(
                np.random.normal(0, 1, (S.shape[0], ST.shape[0])), 
                device=device, requires_grad=True, dtype=torch.float32
            )
            self.mask = None
        else:
            self.cluster_matrix = torch.tensor(
                cluster_matrix, device=device, requires_grad=True, dtype=torch.float32
            )
            self.mask = self.cluster_matrix.clone().detach() < 1
        
        self.device = device
        self.x_length = int(ST.obs['x'].max() - ST.obs['x'].min() + 1)
        self.y_length = int(ST.obs['y'].max() - ST.obs['y'].min() + 1)
        self.x_index = torch.tensor(
            ST.obs['x'] - ST.obs['x'].min(), device=device, dtype=torch.int32
        )
        self.y_index = torch.tensor(
            ST.obs['y'] - ST.obs['y'].min(), device=device, dtype=torch.int32
        )

    def _generate_Xstrips(self, cluster_probs):
        """
        Generate X matrix based on cluster probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: X matrix.
        """
        mapping_x = torch.zeros((self.x_length, self.S.shape[0]), device=self.device)
        mapping_x.index_add_(0, self.x_index, cluster_probs.T)
        return mapping_x

    def _generate_Ystrips(self, cluster_probs):
        """
        Generate Y matrix based on cluster probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: Y matrix.
        """
        mapping_y = torch.zeros((self.y_length, self.S.shape[0]), device=self.device)
        mapping_y.index_add_(0, self.y_index, cluster_probs.T)
        return mapping_y
    
    def _loss_fn(self, verbose=True):
        """
        Computes the loss function for the cell-to-spot mapping.
        
        Args:
            verbose (bool): If True, prints the loss.
        
        Returns:
            tuple: Total loss, X-axis main loss terms, Y-axis main loss terms, and mask sum.
        """
        cluster_probs = softmax(self.cluster_matrix, dim=1)
        mask_sum = torch.masked_select(cluster_probs, self.mask).mean() if self.mask is not None else 0
        
        Mx = self._generate_Xstrips(cluster_probs)
        My = self._generate_Ystrips(cluster_probs)
        
        Gx_pred = torch.matmul(Mx, self.S)
        Gy_pred = torch.matmul(My, self.S)
        
        gx1 = self.lambda_gx1 * cosine_similarity(Gx_pred, self.Gx, dim=0).mean()
        gx2 = self.lambda_gx2 * cosine_similarity(Gx_pred, self.Gx, dim=1).mean()
        gy1 = self.lambda_gy1 * cosine_similarity(Gy_pred, self.Gy, dim=0).mean()
        gy2 = self.lambda_gy2 * cosine_similarity(Gy_pred, self.Gy, dim=1).mean()
        
        cos_x = gx1 + gx2
        cos_y = gy1 + gy2
        
        total_loss = -cos_x - cos_y + mask_sum
        
        if verbose:
            print(f"Total loss: {total_loss.item():.3f}")
        
        return total_loss, gx1, gx2, gy1, gy2, mask_sum

    def fit(self, num_epochs=500,learning_rate=0.1, print_each=100):
        """
        Runs the optimization process and returns the mapping matrix.
        
        Args:
            num_epochs (int): Number of epochs to train.
            learning_rate (float): Learning rate for the optimizer.
            print_each (int): Frequency of loss printing.
        
        Returns:
            ndarray: Optimized mapping matrix (number_cells x number_spots).
        """
        torch.manual_seed(seed=self.random_state)
        optimizer = torch.optim.Adam([self.cluster_matrix], lr=learning_rate)
        
        for t in range(num_epochs):
            loss, gx1, gx2, gy1, gy2, mask_sum = self._loss_fn(verbose=(print_each is not None and t % print_each == 0))
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        
        with torch.no_grad():
            spot_matrix = softmax(self.cluster_matrix, dim=1).cpu().numpy()
        
        if self.device == 'cuda':
            torch.cuda.empty_cache()
        return spot_matrix

class Cell2VisiumSpots:
    """
    Maps each cell to its visium spot based on multi-layer spatial constraints.
    
    Attributes:
        S (torch.Tensor): Single-cell expression matrix.
        ST (AnnData): Spatial data with 3D coordinates.
        Gx (torch.Tensor): X-axis spatial expression matrix.
        Gy (torch.Tensor): Y-axis spatial expression matrix.
        Gz (torch.Tensor): Z-axis spatial expression matrix.
        lambda_gx1 (float): Weight for the first term in X-axis loss function.
        lambda_gx2 (float): Weight for the second term in X-axis loss function.
        lambda_gy1 (float): Weight for the first term in Y-axis loss function.
        lambda_gy2 (float): Weight for the second term in Y-axis loss function.
        lambda_gz1 (float): Weight for the first term in Z-axis loss function.
        lambda_gz2 (float): Weight for the second term in Z-axis loss function.
        cluster_matrix (torch.Tensor): Initial cluster matrix for optimization.
        device (str): Device to use ('cpu' or 'cuda').
        random_state (int): Random seed for initialization.
    """
    
    def __init__(
        self,
        S,
        ST,
        Gx,
        Gy,
        Gz,
        lambda_gx1=1,
        lambda_gx2=0,
        lambda_gy1=1,
        lambda_gy2=0,
        lambda_gz1=1,
        lambda_gz2=0,
        cluster_matrix=None,
        device="cpu",
        random_state=1597,
    ):
        """
        Initializes the Cell2Spots3D instance.
        
        Args:
            S (ndarray): Single-cell expression matrix.
            ST (AnnData): Spatial data with 3D coordinates.
            Gx (ndarray): X-axis spatial expression matrix.
            Gy (ndarray): Y-axis spatial expression matrix.
            Gz (ndarray): Z-axis spatial expression matrix.
            lambda_gx1 (float): Weight for the first term in X-axis loss function.
            lambda_gx2 (float): Weight for the second term in X-axis loss function.
            lambda_gy1 (float): Weight for the first term in Y-axis loss function.
            lambda_gy2 (float): Weight for the second term in Y-axis loss function.
            lambda_gz1 (float): Weight for the first term in Z-axis loss function.
            lambda_gz2 (float): Weight for the second term in Z-axis loss function.
            cluster_matrix (ndarray or None): Initial cluster matrix for optimization.
            device (str): Device to use ('cpu' or 'cuda').
            random_state (int): Random seed for initialization.
        """
        self.S = torch.tensor(S, device=device, dtype=torch.float32)
        self.ST = ST
        self.Gx = torch.tensor(Gx, device=device, dtype=torch.float32)
        self.Gy = torch.tensor(Gy, device=device, dtype=torch.float32)
        self.Gz = torch.tensor(Gz, device=device, dtype=torch.float32)
        self.lambda_gx1 = lambda_gx1
        self.lambda_gx2 = lambda_gx2
        self.lambda_gy1 = lambda_gy1
        self.lambda_gy2 = lambda_gy2
        self.lambda_gz1 = lambda_gz1
        self.lambda_gz2 = lambda_gz2
        self.random_state = random_state
        np.random.seed(seed=self.random_state)
        
        if cluster_matrix is None:
            self.cluster_matrix = np.random.normal(0, 1, (S.shape[0], ST.shape[0]))
        else:
            self.cluster_matrix = cluster_matrix
        
        self.cluster_matrix = torch.tensor(
            self.cluster_matrix, device=device, requires_grad=True, dtype=torch.float32
        )
        
        self.device = device
        self.x_length = int(ST.obs['x'].max() - ST.obs['x'].min() + 1)
        self.y_length = int(ST.obs['y'].max() - ST.obs['y'].min() + 1)
        self.z_length = int(ST.obs['z'].max() - ST.obs['z'].min() + 1)
        
        self.x_index = torch.tensor(ST.obs['x'] - ST.obs['x'].min(), device=device, dtype=torch.int32)
        self.y_index = torch.tensor(ST.obs['y'] - ST.obs['y'].min(), device=device, dtype=torch.int32)
        self.z_index = torch.tensor(ST.obs['z'] - ST.obs['z'].min(), device=device, dtype=torch.int32)

    def _generate_Xstrips(self, cluster_probs):
        """
        Generate X matrix based on regional probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: X matrix.
        """
        mapping_x = torch.zeros((self.x_length, self.S.shape[0]), device=self.device)
        mapping_x.index_add_(0, self.x_index, cluster_probs.T)
        return mapping_x

    def _generate_Ystrips(self, cluster_probs):
        """
        Generate Y matrix based on regional probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: Y matrix.
        """
        mapping_y = torch.zeros((self.y_length, self.S.shape[0]), device=self.device)
        mapping_y.index_add_(0, self.y_index, cluster_probs.T)
        return mapping_y
    
    def _generate_Zstrips(self, cluster_probs):
        """
        Generate Z matrix based on regional probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: Z matrix.
        """
        mapping_z = torch.zeros((self.z_length, self.S.shape[0]), device=self.device)
        mapping_z.index_add_(0, self.z_index, cluster_probs.T)
        return mapping_z

    def _loss_fn(self, verbose=True):
        """
        Computes the loss function for the cell-to-spot mapping.
        
        Args:
            verbose (bool): If True, prints the loss.
        
        Returns:
            tuple: Total loss, and individual loss terms for X, Y, and Z axes.
        """
        cluster_probs = softmax(self.cluster_matrix, dim=1)
        
        Mx = self._generate_Xstrips(cluster_probs)
        My = self._generate_Ystrips(cluster_probs)
        Mz = self._generate_Zstrips(cluster_probs)
        
        Gx_pred = torch.matmul(Mx, self.S)
        Gy_pred = torch.matmul(My, self.S)
        Gz_pred = torch.matmul(Mz, self.S)
        
        x_gv_term = self.lambda_gx1 * cosine_similarity(Gx_pred, self.Gx, dim=0).mean()
        x_vg_term = self.lambda_gx2 * cosine_similarity(Gx_pred, self.Gx, dim=1).mean()
        y_gv_term = self.lambda_gy1 * cosine_similarity(Gy_pred, self.Gy, dim=0).mean()
        y_vg_term = self.lambda_gy2 * cosine_similarity(Gy_pred, self.Gy, dim=1).mean()
        z_gv_term = self.lambda_gz1 * cosine_similarity(Gz_pred, self.Gz, dim=0).mean()
        z_vg_term = self.lambda_gz2 * cosine_similarity(Gz_pred, self.Gz, dim=1).mean()        
        
        expression_term_x = x_gv_term + x_vg_term
        expression_term_y = y_gv_term + y_vg_term
        expression_term_z = z_gv_term + z_vg_term
        
        total_loss = -expression_term_x - expression_term_y - expression_term_z
        
        if verbose:
            print(f"Total loss: {total_loss.item():.3f}")
        
        return (
            total_loss,
            x_gv_term,
            x_vg_term,
            y_gv_term,
            y_vg_term,
            z_gv_term,
            z_vg_term
        )

    def fit(self, num_epochs=500,learning_rate=0.1, print_each=100):
        """
        Runs the optimization process and returns the mapping matrix.
        
        Args:
            num_epochs (int): Number of epochs to train.
            learning_rate (float): Learning rate for the optimizer.
            print_each (int): Frequency of loss printing.
        
        Returns:
            ndarray: Optimized mapping matrix (number_cells x number_spots).
        """
        torch.manual_seed(seed=self.random_state)
        optimizer = torch.optim.Adam([self.cluster_matrix], lr=learning_rate)
        
        for t in range(num_epochs):
            loss, *loss_terms = self._loss_fn(verbose=(print_each is not None and t % print_each == 0))
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        
        with torch.no_grad():
            spot_matrix = softmax(self.cluster_matrix, dim=1).cpu().numpy()
        
        if self.device == 'cuda':
            torch.cuda.empty_cache()
        return spot_matrix


class Cell2SpotsStrips:
    """
    Maps each cell to its spot using strips bulk RNA-seq data.
    
    Attributes:
        S (torch.Tensor): Single-cell expression matrix.
        Gx (torch.Tensor): X-axis spatial expression matrix.
        Gy (torch.Tensor): Y-axis spatial expression matrix.
        lambda_gx1 (float): Weight for the first term in X-axis loss function.
        lambda_gx2 (float): Weight for the second term in X-axis loss function.
        lambda_gy1 (float): Weight for the first term in Y-axis loss function.
        lambda_gy2 (float): Weight for the second term in Y-axis loss function.
        cluster_matrix (torch.Tensor): Initial cluster matrix for optimization.
        device (str): Device to use ('cpu' or 'cuda').
        random_state (int): Random seed for initialization.
        mask (torch.Tensor): Mask for filtering out low values.
        x_length (int): Number of X-axis spots.
        y_length (int): Number of Y-axis spots.
        x_index (torch.Tensor): Indices for X-axis spots.
        y_index (torch.Tensor): Indices for Y-axis spots.
    """
    
    def __init__(
        self,
        S,
        Gx,
        Gy,
        lambda_gx1=1,
        lambda_gx2=0,
        lambda_gy1=1,
        lambda_gy2=0,
        cluster_matrix=None,
        device="cpu",
        random_state=1597,
    ):
        """
        Initializes the Cell2SpotsStrips instance.
        
        Args:
            S (ndarray): Single-cell expression matrix.
            Gx (ndarray): X-axis spatial expression matrix.
            Gy (ndarray): Y-axis spatial expression matrix.
            lambda_gx1 (float): Weight for the first term in X-axis loss function.
            lambda_gx2 (float): Weight for the second term in X-axis loss function.
            lambda_gy1 (float): Weight for the first term in Y-axis loss function.
            lambda_gy2 (float): Weight for the second term in Y-axis loss function.
            cluster_matrix (ndarray or None): Initial cluster matrix for optimization.
            device (str): Device to use ('cpu' or 'cuda').
            random_state (int): Random seed for initialization.
        """
        self.S = torch.tensor(S, device=device, dtype=torch.float32)
        self.Gx = torch.tensor(Gx, device=device, dtype=torch.float32)
        self.Gy = torch.tensor(Gy, device=device, dtype=torch.float32)
        self.lambda_gx1 = lambda_gx1
        self.lambda_gx2 = lambda_gx2
        self.lambda_gy1 = lambda_gy1
        self.lambda_gy2 = lambda_gy2
        self.random_state = random_state
        np.random.seed(seed=self.random_state)
        
        if cluster_matrix is None:
            self.cluster_matrix = np.random.normal(0, 1, (S.shape[0], Gx.shape[1]))
        else:
            self.cluster_matrix = cluster_matrix
        
        self.cluster_matrix = torch.tensor(
            self.cluster_matrix, device=device, requires_grad=True, dtype=torch.float32
        )
        
        self.device = device
        self.x_length = Gx.shape[0]
        self.y_length = Gy.shape[0]
        
        # Assuming ST is available with x and y coordinates; otherwise, provide default indices
        self.x_index = torch.arange(self.x_length, device=device, dtype=torch.int32)
        self.y_index = torch.arange(self.y_length, device=device, dtype=torch.int32)
        
        # Mask initialization to filter out low values
        self.mask = self.cluster_matrix.clone().detach() < 1

    def _generate_Xstrips(self, cluster_probs):
        """
        Generate X matrix based on cluster probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: X matrix.
        """
        mapping_x = torch.zeros((self.x_length, self.S.shape[0]), device=self.device)
        mapping_x.index_add_(0, self.x_index, cluster_probs.T)
        return mapping_x

    def _generate_Ystrips(self, cluster_probs):
        """
        Generate Y matrix based on cluster probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: Y matrix.
        """
        mapping_y = torch.zeros((self.y_length, self.S.shape[0]), device=self.device)
        mapping_y.index_add_(0, self.y_index, cluster_probs.T)
        return mapping_y
    
    def _loss_fn(self, verbose=True):
        """
        Computes the loss function for cell-to-spot mapping.
        
        Args:
            verbose (bool): If True, prints the loss.
        
        Returns:
            tuple: Total loss, and individual loss terms for X and Y axes, and mask sum.
        """
        cluster_probs = softmax(self.cluster_matrix, dim=1)
        
        mask_sum = torch.masked_select(cluster_probs, self.mask).mean()
        Mx = self._generate_Xstrips(cluster_probs)
        My = self._generate_Ystrips(cluster_probs)
        
        dx_pred = torch.log(
            Mx.T.sum(axis=0) / Mx.shape[1]
        ) 
        dy_pred = torch.log(
            My.T.sum(axis=0) / My.shape[1]
        )
        
        Gx_pred = torch.matmul(Mx, self.S)
        Gy_pred = torch.matmul(My, self.S)
        
        gx1 = self.lambda_gx1 * cosine_similarity(Gx_pred, self.Gx, dim=0).mean()
        gx2 = self.lambda_gx2 * cosine_similarity(Gx_pred, self.Gx, dim=1).mean()
        gy1 = self.lambda_gy1 * cosine_similarity(Gy_pred, self.Gy, dim=0).mean()
        gy2 = self.lambda_gy2 * cosine_similarity(Gy_pred, self.Gy, dim=1).mean()
        
        cos_x = gx1 + gx2
        cos_y = gy1 + gy2
        
        main_x_loss1 = (gx1 / self.lambda_gx1).tolist()
        main_y_loss1 = (gy1 / self.lambda_gy1).tolist()
        main_x_loss2 = (gx2 / self.lambda_gx2).tolist()
        main_y_loss2 = (gy2 / self.lambda_gy2).tolist()
        
        total_loss = -cos_x - cos_y + mask_sum
        
        if verbose:
            print(f"Total loss: {total_loss.item():.3f}")
        
        return (
            total_loss,
            main_x_loss1,
            main_y_loss1,
            mask_sum
        )

    def fit(self, num_epochs=500,learning_rate=0.1, print_each=100):
        """
        Runs the optimization process and returns the mapping matrix.
        
        Args:
            num_epochs (int): Number of epochs to train.
            learning_rate (float): Learning rate for the optimizer.
            print_each (int): Frequency of loss printing.
        
        Returns:
            ndarray: Optimized mapping matrix (number_cells x number_spots).
        """
        torch.manual_seed(seed=self.random_state)
        optimizer = torch.optim.Adam([self.cluster_matrix], lr=learning_rate)
        
        for t in range(num_epochs):
            if print_each is None or t % print_each != 0:
                run_loss = self._loss_fn(verbose=False)
            else:
                run_loss = self._loss_fn(verbose=True)
            
            loss = run_loss[0]
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        
        with torch.no_grad():
            spot_matrix = softmax(self.cluster_matrix, dim=1).cpu().numpy()
        
        if self.device == 'cuda':
            torch.cuda.empty_cache()
        
        return spot_matrix


class Cell2SpotsEnhance:
    """
    Maps each cell to its spot using multi-layer custom regional information.
    
    Attributes:
        S (torch.Tensor): Single-cell expression matrix.
        ST (AnnData): Spatial transcriptomics data.
        Gx (torch.Tensor): X-axis spatial expression matrix.
        Gy (torch.Tensor): Y-axis spatial expression matrix.
        lambda_gx1 (float): Weight for the first term in X-axis loss function.
        lambda_gx2 (float): Weight for the second term in X-axis loss function.
        lambda_gy1 (float): Weight for the first term in Y-axis loss function.
        lambda_gy2 (float): Weight for the second term in Y-axis loss function.
        cluster_matrix (torch.Tensor): Initial cluster matrix for optimization.
        custom_regions (list): List of custom region labels for regional loss computation.
        lambda_r (list): List of weights for the custom regions.
        device (str): Device to use ('cpu' or 'cuda').
        random_state (int): Random seed for initialization.
        mask (torch.Tensor): Mask for filtering out low values.
        x_length (int): Number of X-axis spots.
        y_length (int): Number of Y-axis spots.
        x_index (torch.Tensor): Indices for X-axis spots.
        y_index (torch.Tensor): Indices for Y-axis spots.
    """
    
    def __init__(
        self,
        S,
        ST,
        Gx,
        Gy,
        lambda_gx1=1,
        lambda_gx2=0,
        lambda_gy1=1,
        lambda_gy2=0,
        cluster_matrix=None,
        custom_regions=None,
        lambda_r=None,
        device="cpu",
        random_state=1597,
    ):
        """
        Initializes the Cell2SpotsEnhance instance.
        
        Args:
            S (ndarray): Single-cell expression matrix.
            ST (AnnData): Spatial transcriptomics data.
            Gx (ndarray): X-axis spatial expression matrix.
            Gy (ndarray): Y-axis spatial expression matrix.
            lambda_gx1 (float): Weight for the first term in X-axis loss function.
            lambda_gx2 (float): Weight for the second term in X-axis loss function.
            lambda_gy1 (float): Weight for the first term in Y-axis loss function.
            lambda_gy2 (float): Weight for the second term in Y-axis loss function.
            cluster_matrix (ndarray or None): Initial cluster matrix for optimization.
            custom_regions (list or None): Custom regions for regional loss computation.
            lambda_r (list or None): Weights for custom regions.
            device (str): Device to use ('cpu' or 'cuda').
            random_state (int): Random seed for initialization.
        """
        self.S = torch.tensor(S, device=device, dtype=torch.float32)
        self.ST = ST
        self.Gx = torch.tensor(Gx, device=device, dtype=torch.float32)
        self.Gy = torch.tensor(Gy, device=device, dtype=torch.float32)
        self.lambda_gx1 = lambda_gx1
        self.lambda_gx2 = lambda_gx2
        self.lambda_gy1 = lambda_gy1
        self.lambda_gy2 = lambda_gy2
        self.random_state = random_state
        self.lambda_r = lambda_r if lambda_r is not None else [1] * (len(custom_regions) if custom_regions else 0)
        np.random.seed(seed=self.random_state)
        
        if cluster_matrix is None:
            self.cluster_matrix = np.random.normal(0, 1, (S.shape[0], ST.shape[0]))
            self.cluster_matrix = torch.tensor(
                self.cluster_matrix, device=device, requires_grad=True, dtype=torch.float32
            )
            self.mask = None
        else:
            self.cluster_matrix = cluster_matrix
            self.cluster_matrix = torch.tensor(
                self.cluster_matrix, device=device, requires_grad=True, dtype=torch.float32
            )
            # self.rawM = self.M.clone().detach() > 0
            self.mask = self.cluster_matrix.clone().detach() < 1
        if custom_regions is None:
            self.region = None
        else:
            self.region = custom_regions
        
        self.device = device
        self.x_length = np.int32(self.ST.obs['x'].max() - self.ST.obs['x'].min() + 1)
        self.y_length = np.int32(self.ST.obs['y'].max() - self.ST.obs['y'].min() + 1)
        self.x_index = self.ST.obs['x'] - self.ST.obs['x'].min()
        self.x_index = torch.tensor(
            self.x_index, device=device, requires_grad=False, dtype=torch.int32
        )
        self.y_index = self.ST.obs['y'] - self.ST.obs['y'].min()
        self.y_index = torch.tensor(
            self.y_index, device=device, requires_grad=False, dtype=torch.int32
        )
        
        self.region = custom_regions

    def _generate_region(self, cluster_probs, region_labels):
        """
        Generate region matrices for custom regions.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
            region_labels (list): List of custom region labels.
        
        Returns:
            tuple: List of region matrices and list of region expression matrices.
        """
        mapping_region = []
        region_exp = []
        ST_exp = torch.tensor(self.ST.X, device=self.device, dtype=torch.float32)
        
        for label in region_labels:
            st_region = torch.tensor(self.ST.obs[label].values, device=self.device, dtype=torch.float32)
            mapping_region.append(st_region.T @ cluster_probs.T)
            region_exp.append(st_region.T @ ST_exp)
        
        return mapping_region, region_exp

    def _generate_Xstrips(self, cluster_probs):
        """
        Generate X matrix based on cluster probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: X matrix.
        """
        mapping_x = torch.zeros((self.x_length, self.S.shape[0]), device=self.device)
        mapping_x.index_add_(0, self.x_index, cluster_probs.T)
        return mapping_x

    def _generate_Ystrips(self, cluster_probs):
        """
        Generate Y matrix based on cluster probabilities.
        
        Args:
            cluster_probs (torch.Tensor): Cluster probabilities.
        
        Returns:
            torch.Tensor: Y matrix.
        """
        mapping_y = torch.zeros((self.y_length, self.S.shape[0]), device=self.device)
        mapping_y.index_add_(0, self.y_index, cluster_probs.T)
        return mapping_y

    def _loss_fn(self, verbose=True):
        """
        Computes the loss function for cell-to-spot mapping.
        
        Args:
            verbose (bool): If True, prints the loss.
        
        Returns:
            tuple: Total loss and individual loss terms for X and Y axes, and mask sum.
        """
        cluster_probs = softmax(self.cluster_matrix, dim=1)
        
        mask_sum = torch.masked_select(cluster_probs, self.mask).mean() if self.mask is not None else 0
        
        Mx = self._generate_Xstrips(cluster_probs)
        My = self._generate_Ystrips(cluster_probs)

        dx_pred = torch.log(Mx.T.sum(axis=0) / Mx.shape[1])
        dy_pred = torch.log(My.T.sum(axis=0) / My.shape[1])
        
        Gx_pred = torch.matmul(Mx, self.S)
        Gy_pred = torch.matmul(My, self.S)
        
        gx1 = self.lambda_gx1 * cosine_similarity(Gx_pred, self.Gx, dim=0).mean()
        gx2 = self.lambda_gx2 * cosine_similarity(Gx_pred, self.Gx, dim=1).mean()
        gy1 = self.lambda_gy1 * cosine_similarity(Gy_pred, self.Gy, dim=0).mean()
        gy2 = self.lambda_gy2 * cosine_similarity(Gy_pred, self.Gy, dim=1).mean()
        
        cos_x = gx1 + gx2
        cos_y = gy1 + gy2
        
        main_x_loss1 = (gx1 / self.lambda_gx1).tolist()
        main_y_loss1 = (gy1 / self.lambda_gy1).tolist()
        main_x_loss2 = (gx2 / self.lambda_gx2).tolist()
        main_y_loss2 = (gy2 / self.lambda_gy2).tolist()
        
        total_loss = -cos_x - cos_y + mask_sum
        
        if self.region:
            Mr, Er = self._generate_region(cluster_probs, self.region)
            for i in range(len(self.region)):
                Gr_pred = torch.matmul(Mr[i], self.S)
                region_loss = self.lambda_r[i] * cosine_similarity(Gr_pred, Er[i], dim=0).mean()
                total_loss -= region_loss

        if verbose:
            print(f"Total loss: {total_loss.item():.3f}")
        
        return (
            total_loss,
            main_x_loss1,
            main_y_loss1,
            mask_sum
        )

    def fit(self,  num_epochs, learning_rate=0.1, print_each=100):
        """
        Runs the optimization process and returns the mapping matrix.
        
        Args:
            num_epochs (int): Number of epochs to train.
            learning_rate (float): Learning rate for the optimizer.
            print_each (int): Frequency of loss printing.
        
        Returns:
            ndarray: Optimized mapping matrix (number_cells x number_spots).
        """
        torch.manual_seed(seed=self.random_state)
        optimizer = torch.optim.Adam([self.cluster_matrix], lr=learning_rate)
        
        for t in range(num_epochs):
            if print_each is None or t % print_each != 0:
                run_loss = self._loss_fn(verbose=False)
            else:
                run_loss = self._loss_fn(verbose=True)
                
            loss = run_loss[0]
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        
        with torch.no_grad():
            spot_matrix = softmax(self.cluster_matrix, dim=1).cpu().numpy()
            torch.cuda.empty_cache()
        
        return spot_matrix
