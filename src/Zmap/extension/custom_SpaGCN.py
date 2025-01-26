import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import issparse
import torch
from sklearn.decomposition import PCA
import math
import torch.nn as nn
import torch.nn.functional as F
from torch.nn.parameter import Parameter
from torch.nn.modules.module import Module
from sklearn.cluster import KMeans
import torch.optim as optim
from random import shuffle
import scipy
from anndata import AnnData,read_csv,read_text,read_mtx
import random
import numba

@numba.njit("f4(f4[:], f4[:])")
def euclid_dist(t1,t2):
	sum=0
	for i in range(t1.shape[0]):
		sum+=(t1[i]-t2[i])**2
	return np.sqrt(sum)

@numba.njit("f4[:,:](f4[:,:])", parallel=True, nogil=True)
def pairwise_distance(X):
	n=X.shape[0]
	adj=np.empty((n, n), dtype=np.float32)
	for i in numba.prange(n):
		for j in numba.prange(n):
			adj[i][j]=euclid_dist(X[i], X[j])
	return adj

def extract_color(x_pixel=None, y_pixel=None, image=None, beta=49):
	#beta to control the range of neighbourhood when calculate grey vale for one spot
	beta_half=round(beta/2)
	g=[]
	for i in range(len(x_pixel)):
		max_x=image.shape[0]
		max_y=image.shape[1]
		nbs=image[max(0,x_pixel[i]-beta_half):min(max_x,x_pixel[i]+beta_half+1),max(0,y_pixel[i]-beta_half):min(max_y,y_pixel[i]+beta_half+1)]
		g.append(np.mean(np.mean(nbs,axis=0),axis=0))
	c0, c1, c2=[], [], []
	for i in g:
		c0.append(i[0])
		c1.append(i[1])
		c2.append(i[2])
	c0=np.array(c0)
	c1=np.array(c1)
	c2=np.array(c2)
	c3=(c0*np.var(c0)+c1*np.var(c1)+c2*np.var(c2))/(np.var(c0)+np.var(c1)+np.var(c2))
	return c3

def calculate_adj_matrix(x, y, x_pixel=None, y_pixel=None, image=None, beta=49, alpha=1, histology=True):
	#x,y,x_pixel, y_pixel are lists
	if histology:
		assert (x_pixel is not None) & (x_pixel is not None) & (image is not None)
		assert (len(x)==len(x_pixel)) & (len(y)==len(y_pixel))
		# print("Calculateing adj matrix using histology image...")
		#beta to control the range of neighbourhood when calculate grey vale for one spot
		#alpha to control the color scale
		beta_half=round(beta/2)
		g=[]
		for i in range(len(x_pixel)):
			max_x=image.shape[0]
			max_y=image.shape[1]
			nbs=image[max(0,x_pixel[i]-beta_half):min(max_x,x_pixel[i]+beta_half+1),max(0,y_pixel[i]-beta_half):min(max_y,y_pixel[i]+beta_half+1)]
			g.append(np.mean(np.mean(nbs,axis=0),axis=0))
		c0, c1, c2=[], [], []
		for i in g:
			c0.append(i[0])
			c1.append(i[1])
			c2.append(i[2])
		c0=np.array(c0)
		c1=np.array(c1)
		c2=np.array(c2)
		print("Var of c0,c1,c2 = ", np.var(c0),np.var(c1),np.var(c2))
		c3=(c0*np.var(c0)+c1*np.var(c1)+c2*np.var(c2))/(np.var(c0)+np.var(c1)+np.var(c2))
		c4=(c3-np.mean(c3))/np.std(c3)
		z_scale=np.max([np.std(x), np.std(y)])*alpha
		z=c4*z_scale
		z=z.tolist()
		print("Var of x,y,z = ", np.var(x),np.var(y),np.var(z))
		X=np.array([x, y, z]).T.astype(np.float32)
	else:
		# print("Calculateing adj matrix using xy only...")
		X=np.array([x, y]).T.astype(np.float32)
	return pairwise_distance(X)
def prefilter_cells(adata,min_counts=None,max_counts=None,min_genes=200,max_genes=None):
    if min_genes is None and min_counts is None and max_genes is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp=np.asarray([True]*adata.shape[0],dtype=bool)
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_cells(adata.X,min_genes=min_genes)[0]) if min_genes is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_cells(adata.X,max_genes=max_genes)[0]) if max_genes is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_cells(adata.X,min_counts=min_counts)[0]) if min_counts is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_cells(adata.X,max_counts=max_counts)[0]) if max_counts is not None  else id_tmp
    adata._inplace_subset_obs(id_tmp)
    adata.raw=sc.pp.log1p(adata,copy=True) #check the rowname 
    # print("the var_names of adata.raw: adata.raw.var_names.is_unique=:",adata.raw.var_names.is_unique)
   

def prefilter_genes(adata,min_counts=None,max_counts=None,min_cells=10,max_cells=None):
    if min_cells is None and min_counts is None and max_cells is None and max_counts is None:
        raise ValueError('Provide one of min_counts, min_genes, max_counts or max_genes.')
    id_tmp=np.asarray([True]*adata.shape[1],dtype=bool)
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_cells=min_cells)[0]) if min_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_cells=max_cells)[0]) if max_cells is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,min_counts=min_counts)[0]) if min_counts is not None  else id_tmp
    id_tmp=np.logical_and(id_tmp,sc.pp.filter_genes(adata.X,max_counts=max_counts)[0]) if max_counts is not None  else id_tmp
    adata._inplace_subset_var(id_tmp)


def prefilter_specialgenes(adata,Gene1Pattern="ERCC",Gene2Pattern="MT-"):
    id_tmp1=np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names],dtype=bool)
    id_tmp2=np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names],dtype=bool)
    id_tmp=np.logical_and(id_tmp1,id_tmp2)
    adata._inplace_subset_var(id_tmp)


def calculate_p(adj, l):
    adj_exp=np.exp(-1*(adj**2)/(2*(l**2)))
    return np.mean(np.sum(adj_exp,1))-1

def test_l(adj, list_l):
    for l in list_l:
        print("l is ",str(l),"Percentage of total expression contributed by neighborhoods:",calculate_p(adj, l))

def find_l(p, adj, start=0.5, end=2,sep=0.01, tol=0.01):
    for l in np.arange(start, end, sep):
        q=calculate_p(adj, l)
        print("L=", str(l), "P=", str(round(q,5)))
        if np.abs(p-q)<=tol:
            return l
    # print("l not found, try bigger range or smaller sep!")

def search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100):
    run=0
    p_low=calculate_p(adj, start)
    p_high=calculate_p(adj, end)
    if p_low>p+tol:
        # print("l not found, try smaller start point.")
        return None
    elif p_high<p-tol:
        # print("l not found, try bigger end point.")
        return None
    elif  np.abs(p_low-p) <=tol:
        # print("recommended l = ", str(start))
        return start
    elif  np.abs(p_high-p) <=tol:
        # print("recommended l = ", str(end))
        return end
    while (p_low+tol)<p<(p_high-tol):
        run+=1
        # print("Run "+str(run)+": l ["+str(start)+", "+str(end)+"], p ["+str(p_low)+", "+str(p_high)+"]")
        if run >max_run:
            # print("Exact l not found, closest values are:\n"+"l="+str(start)+": "+"p="+str(p_low)+"\nl="+str(end)+": "+"p="+str(p_high))
            return None
        mid=(start+end)/2
        p_mid=calculate_p(adj, mid)
        if np.abs(p_mid-p)<=tol:
            # print("recommended l = ", str(mid))
            return mid
        if p_mid<=p:
            start=mid
            p_low=p_mid
        else:
            end=mid
            p_high=p_mid

def count_nbr(target_cluster,cell_id, x, y, pred, radius):
    adj_2d=calculate_adj_matrix(x=x,y=y, histology=False)
    cluster_num = dict()
    df = {'cell_id': cell_id, 'x': x, "y":y, "pred":pred}
    df = pd.DataFrame(data=df)
    df.index=df['cell_id']
    target_df=df[df["pred"]==target_cluster]
    row_index=0
    num_nbr=[]
    for index, row in target_df.iterrows():
        x=row["x"]
        y=row["y"]
        tmp_nbr=df[((df["x"]-x)**2+(df["y"]-y)**2)<=(radius**2)]
        num_nbr.append(tmp_nbr.shape[0])
    return np.mean(num_nbr)

def search_radius(target_cluster,cell_id, x, y, pred, start, end, num_min=8, num_max=15,  max_run=100):
    run=0
    num_low=count_nbr(target_cluster,cell_id, x, y, pred, start)
    num_high=count_nbr(target_cluster,cell_id, x, y, pred, end)
    if num_min<=num_low<=num_max:
        # print("recommended radius = ", str(start))
        return start
    elif num_min<=num_high<=num_max:
        # print("recommended radius = ", str(end))
        return end
    elif num_low>num_max:
        # print("Try smaller start.")
        return None
    elif num_high<num_min:
        # print("Try bigger end.")
        return None
    while (num_low<num_min) and (num_high>num_min):
        run+=1
        # print("Run "+str(run)+": radius ["+str(start)+", "+str(end)+"], num_nbr ["+str(num_low)+", "+str(num_high)+"]")
        if run >max_run:
            # print("Exact radius not found, closest values are:\n"+"radius="+str(start)+": "+"num_nbr="+str(num_low)+"\nradius="+str(end)+": "+"num_nbr="+str(num_high))
            return None
        mid=(start+end)/2
        num_mid=count_nbr(target_cluster,cell_id, x, y, pred, mid)
        if num_min<=num_mid<=num_max:
            # print("recommended radius = ", str(mid), "num_nbr="+str(num_mid))
            return mid
        if num_mid<num_min:
            start=mid
            num_low=num_mid
        elif num_mid>num_max:
            end=mid
            num_high=num_mid

def find_neighbor_clusters(target_cluster,cell_id, x, y, pred,radius, ratio=1/2):
    cluster_num = dict()
    for i in pred:
        cluster_num[i] = cluster_num.get(i, 0) + 1
    df = {'cell_id': cell_id, 'x': x, "y":y, "pred":pred}
    df = pd.DataFrame(data=df)
    df.index=df['cell_id']
    target_df=df[df["pred"]==target_cluster]
    nbr_num={}
    row_index=0
    num_nbr=[]
    for index, row in target_df.iterrows():
        x=row["x"]
        y=row["y"]
        tmp_nbr=df[((df["x"]-x)**2+(df["y"]-y)**2)<=(radius**2)]
        #tmp_nbr=df[(df["x"]<x+radius) & (df["x"]>x-radius) & (df["y"]<y+radius) & (df["y"]>y-radius)]
        num_nbr.append(tmp_nbr.shape[0])
        for p in tmp_nbr["pred"]:
            nbr_num[p]=nbr_num.get(p,0)+1
    del nbr_num[target_cluster]
    nbr_num_back=nbr_num.copy() #Backup
    nbr_num=[(k, v)  for k, v in nbr_num.items() if v>(ratio*cluster_num[k])]
    nbr_num.sort(key=lambda x: -x[1])
    # print("radius=", radius, "average number of neighbors for each spot is", np.mean(num_nbr))
    # print(" Cluster",target_cluster, "has neighbors:")
    for t in nbr_num:
        print("Dmain ", t[0], ": ",t[1])
    ret=[t[0] for t in nbr_num]
    if len(ret)==0:
        nbr_num_back=[(k, v)  for k, v in nbr_num_back.items()]
        nbr_num_back.sort(key=lambda x: -x[1])
        ret=[nbr_num_back[0][0]]
        # print("No neighbor domain found, only return one potential neighbor domain:",ret)
        # print("Try bigger radius or smaller ratio.")
    return ret


def rank_genes_groups(input_adata, target_cluster,nbr_list, label_col, adj_nbr=True, log=False):
    if adj_nbr:
        nbr_list=nbr_list+[target_cluster]
        adata=input_adata[input_adata.obs[label_col].isin(nbr_list)]
    else:
        adata=input_adata.copy()
    adata.var_names_make_unique()
    adata.obs["target"]=((adata.obs[label_col]==target_cluster)*1).astype('category')
    sc.tl.rank_genes_groups(adata, groupby="target",reference="rest", n_genes=adata.shape[1],method='wilcoxon')
    pvals_adj=[i[0] for i in adata.uns['rank_genes_groups']["pvals_adj"]]
    genes=[i[1] for i in adata.uns['rank_genes_groups']["names"]]
    if issparse(adata.X):
        obs_tidy=pd.DataFrame(adata.X.A)
    else:
        obs_tidy=pd.DataFrame(adata.X)
    obs_tidy.index=adata.obs["target"].tolist()
    obs_tidy.columns=adata.var.index.tolist()
    obs_tidy=obs_tidy.loc[:,genes]
    # 1. compute mean value
    mean_obs = obs_tidy.groupby(level=0).mean()
    # 2. compute fraction of cells having value >0
    obs_bool = obs_tidy.astype(bool)
    fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()
    # compute fold change.
    if log: #The adata already logged
        fold_change=np.exp((mean_obs.loc[1] - mean_obs.loc[0]).values)
    else:
        fold_change = (mean_obs.loc[1] / (mean_obs.loc[0]+ 1e-9)).values
    df = {'genes': genes, 'in_group_fraction': fraction_obs.loc[1].tolist(), "out_group_fraction":fraction_obs.loc[0].tolist(),"in_out_group_ratio":(fraction_obs.loc[1]/fraction_obs.loc[0]).tolist(),"in_group_mean_exp": mean_obs.loc[1].tolist(), "out_group_mean_exp": mean_obs.loc[0].tolist(),"fold_change":fold_change.tolist(), "pvals_adj":pvals_adj}
    df = pd.DataFrame(data=df)
    return df

def relative_func(expres):
    #expres: an array counts expression for a gene
    maxd = np.max(expres) - np.min(expres)
    min_exp=np.min(expres)
    rexpr = (expres - min_exp)/maxd
    return rexpr

def plot_relative_exp(input_adata, gene, x_name, y_name,color,use_raw=False, spot_size=200000):
    adata=input_adata.copy()
    if use_raw:
        X=adata.raw.X
    else:
        X=adata.X
    if issparse(X):
        X=pd.DataFrame(X.A)
    else:
        X=pd.DataFrame(X)
    X.index=adata.obs.index
    X.columns=adata.var.index
    rexpr=relative_func(X.loc[:,gene])
    adata.obs["rexpr"]=rexpr
    fig=sc.pl.scatter(adata,x=x_name,y=y_name,color="rexpr",title=gene+"_rexpr",color_map=color,show=False,size=spot_size/adata.shape[0])
    return fig

def plot_log_exp(input_adata, gene, x_name, y_name,color,use_raw=False):
    adata=input_adata.copy()
    if use_raw:
        X=adata.X
    else:
        X=adata.raw.X
    if issparse(X):
        X=pd.DataFrame(X.A)
    else:
        X=pd.DataFrame(X)
    X.index=adata.obs.index
    X.columns=adata.var.index
    adata.obs["log"]=np.log((X.loc[:,gene]+1).tolist())
    fig=sc.pl.scatter(adata,x=x_name,y=y_name,color="log",title=gene+"_log",color_map=color,show=False,size=200000/adata.shape[0])
    return fig

def detect_subclusters(cell_id, x, y, pred, target_cluster, radius=3, res=0.2):
    df = {'cell_id': cell_id, 'x': x, "y":y, "pred":pred}
    df = pd.DataFrame(data=df)
    df.index=df['cell_id']
    target_df=df[df["pred"]==target_cluster]
    nbr=np.zeros([target_df.shape[0],len(set(df["pred"]))],dtype=int)
    num_nbr=[]
    row_index=0
    for index, row in target_df.iterrows():
        x=row["x"]
        y=row["y"]
        tmp_nbr=df[(df["x"]<x+radius) & (df["x"]>x-radius) & (df["y"]<y+radius) & (df["y"]>y-radius)]
        num_nbr.append(tmp_nbr.shape[0])
        for p in tmp_nbr["pred"]:
            nbr[row_index,int(p)]+=1
        row_index+=1
    #Minus out the cell itself
    nbr[:,target_cluster]=nbr[:,target_cluster]-1
    nbr=sc.AnnData(nbr)
    sc.pp.neighbors(nbr, n_neighbors=10)
    sc.tl.leiden(nbr,resolution=res)
    sub_cluster=nbr.obs['leiden'].astype(int).to_numpy()
    target_df["sub_cluster"]=sub_cluster
    target_df["sub_cluster"]=target_df["sub_cluster"].astype('category')
    tmp=[]
    for j in df.index:
        if j in target_df.index:
            tmp.append(target_df.loc[j,"sub_cluster"])
        else:
            tmp.append("-1")
    #ret = {'cell_id': cell_id, 'sub_cluster_'+str(target_cluster): tmp}
    #ret = pd.DataFrame(data=ret)
    #ret.index=ret['cell_id']
    ret=tmp
    return ret

def find_meta_gene(input_adata,
                    pred,
                    target_domain,
                    start_gene,
                    mean_diff=0,
                    early_stop=True,
                    max_iter=5,
                    use_raw=False):
    meta_name=start_gene
    adata=input_adata.copy()
    adata.obs["meta"]=adata.X[:,adata.var.index==start_gene]
    adata.obs["pred"]=pred
    num_non_target=adata.shape[0]
    for i in range(max_iter):
        #Select cells
        tmp=adata[((adata.obs["meta"]>np.mean(adata.obs[adata.obs["pred"]==target_domain]["meta"]))|(adata.obs["pred"]==target_domain))]
        tmp.obs["target"]=((tmp.obs["pred"]==target_domain)*1).astype('category').copy()
        if (len(set(tmp.obs["target"]))<2) or (np.min(tmp.obs["target"].value_counts().values)<5):
            # print("Meta gene is: ", meta_name)
            return meta_name, adata.obs["meta"].tolist()
        #DE
        sc.tl.rank_genes_groups(tmp, groupby="target",reference="rest", n_genes=1,method='wilcoxon')
        adj_g=tmp.uns['rank_genes_groups']["names"][0][0]
        add_g=tmp.uns['rank_genes_groups']["names"][0][1]
        meta_name_cur=meta_name+"+"+add_g+"-"+adj_g
        # print("Add gene: ", add_g)
        # print("Minus gene: ", adj_g)
        #Meta gene
        adata.obs[add_g]=adata.X[:,adata.var.index==add_g]
        adata.obs[adj_g]=adata.X[:,adata.var.index==adj_g]
        adata.obs["meta_cur"]=(adata.obs["meta"]+adata.obs[add_g]-adata.obs[adj_g])
        adata.obs["meta_cur"]=adata.obs["meta_cur"]-np.min(adata.obs["meta_cur"])
        mean_diff_cur=np.mean(adata.obs["meta_cur"][adata.obs["pred"]==target_domain])-np.mean(adata.obs["meta_cur"][adata.obs["pred"]!=target_domain])
        num_non_target_cur=np.sum(tmp.obs["target"]==0)
        if (early_stop==False) | ((num_non_target>=num_non_target_cur) & (mean_diff<=mean_diff_cur)):
            num_non_target=num_non_target_cur
            mean_diff=mean_diff_cur
            # print("Absolute mean change:", mean_diff)
            # print("Number of non-target spots reduced to:",num_non_target)
        else:
            # print("Stopped!", "Previous Number of non-target spots",num_non_target, num_non_target_cur, mean_diff,mean_diff_cur)
            # print("Previous Number of non-target spots",num_non_target, num_non_target_cur, mean_diff,mean_diff_cur)
            # print("Previous Number of non-target spots",num_non_target)
            # print("Current Number of non-target spots",num_non_target_cur)
            # print("Absolute mean change", mean_diff)
            # print("===========================================================================")
            # print("Meta gene: ", meta_name)
            # print("===========================================================================")
            return meta_name, adata.obs["meta"].tolist()
        meta_name=meta_name_cur
        adata.obs["meta"]=adata.obs["meta_cur"]
        # print("===========================================================================")
        # print("Meta gene is: ", meta_name)
        # print("===========================================================================")
    return meta_name, adata.obs["meta"].tolist()


def search_res(adata, adj,embed, l, target_num, start=0.4, step=0.1, tol=5e-3, lr=0.05, max_epochs=10, r_seed=100, t_seed=100, n_seed=100, max_run=10):
    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)
    res=start
    # print("Start at res = ", res, "step = ", step)
    clf=SpaGCN()
    clf.set_l(l)
    clf.train(adata,adj,embed,init_spa=True,init="leiden",res=res, tol=tol, lr=lr, max_epochs=max_epochs)
    y_pred, _=clf.predict()
    old_num=len(set(y_pred))
    # print("Res = ", res, "Num of clusters = ", old_num)
    run=0
    while old_num!=target_num:
        random.seed(r_seed)
        torch.manual_seed(t_seed)
        np.random.seed(n_seed)
        old_sign=1 if (old_num<target_num) else -1
        clf=SpaGCN()
        clf.set_l(l)
        clf.train(adata,adj,embed,init_spa=True,init="leiden",res=res+step*old_sign, tol=tol, lr=lr, max_epochs=max_epochs)
        y_pred, _=clf.predict()
        new_num=len(set(y_pred))
        # print("Res = ", res+step*old_sign, "Num of clusters = ", new_num)
        if new_num==target_num:
            res=res+step*old_sign
            # print("recommended res = ", str(res))
            return res
        new_sign=1 if (new_num<target_num) else -1
        if new_sign==old_sign:
            res=res+step*old_sign
            # print("Res changed to", res)
            old_num=new_num
        else:
            step=step/2
            # print("Step changed to", step)
        if run >max_run:
            # print("Exact resolution not found")
            # print("Recommended res = ", str(res))
            return res
        run+=1
    # print("recommended res = ", str(res))
    return res


def refine(sample_id, pred, dis, shape="hexagon"):
    refined_pred=[]
    pred=pd.DataFrame({"pred": pred}, index=sample_id)
    dis_df=pd.DataFrame(dis, index=sample_id, columns=sample_id)
    if shape=="hexagon":
        num_nbs=6 
    elif shape=="square":
        num_nbs=4
    else:
        print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
    for i in range(len(sample_id)):
        index=sample_id[i]
        dis_tmp=dis_df.loc[index, :].sort_values()
        nbs=dis_tmp[0:num_nbs+1]
        nbs_pred=pred.loc[nbs.index, "pred"]
        self_pred=pred.loc[index, "pred"]
        v_c=nbs_pred.value_counts()
        if (v_c.loc[self_pred]<num_nbs/2) and (np.max(v_c)>num_nbs/2):
            refined_pred.append(v_c.idxmax())
        else:           
            refined_pred.append(self_pred)
    return refined_pred
class GraphConvolution(Module):
    """
    Simple GCN layer, similar to https://arxiv.org/abs/1609.02907
    """

    def __init__(self, in_features, out_features, bias=True):
        super(GraphConvolution, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = Parameter(torch.FloatTensor(in_features, out_features))
        if bias:
            self.bias = Parameter(torch.FloatTensor(out_features))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        stdv = 1. / math.sqrt(self.weight.size(1))
        self.weight.data.uniform_(-stdv, stdv)
        if self.bias is not None:
            self.bias.data.uniform_(-stdv, stdv)

    def forward(self, input, adj):
        support = torch.mm(input, self.weight)
        output = torch.spmm(adj, support)
        if self.bias is not None:
            return output + self.bias
        else:
            return output

    def __repr__(self):
        return self.__class__.__name__ + ' (' \
               + str(self.in_features) + ' -> ' \
               + str(self.out_features) + ')'
class simple_GC_DEC(nn.Module):
    def __init__(self, nfeat, nhid, alpha=0.2):
        super(simple_GC_DEC, self).__init__()
        self.gc = GraphConvolution(nfeat, nhid)
        self.nhid=nhid
        #self.mu determined by the init method
        self.alpha=alpha

    def forward(self, x, adj):
        x=self.gc(x, adj)
        q = 1.0 / ((1.0 + torch.sum((x.unsqueeze(1) - self.mu)**2, dim=2) / self.alpha) + 1e-8)
        q = q**(self.alpha+1.0)/2.0
        q = q / torch.sum(q, dim=1, keepdim=True)
        return x, q

    def loss_function(self, p, q):
        def kld(target, pred):
            return torch.mean(torch.sum(target*torch.log(target/(pred+1e-6)), dim=1))
        loss = kld(p, q)
        return loss

    def target_distribution(self, q):
        #weight = q ** 2 / q.sum(0)
        #return torch.transpose((torch.transpose(weight,0,1) / weight.sum(1)),0,1)e
        p = q**2 / torch.sum(q, dim=0)
        p = p / torch.sum(p, dim=1, keepdim=True)
        return p

    def fit(self, X,adj,  lr=0.001, max_epochs=5000, update_interval=3, trajectory_interval=50,weight_decay=5e-4,opt="sgd",init="leiden",n_neighbors=10,res=0.4,n_clusters=10,init_spa=True,tol=1e-3):
        self.trajectory=[]
        if opt=="sgd":
            optimizer = optim.SGD(self.parameters(), lr=lr, momentum=0.9)
        elif opt=="admin":
            optimizer = optim.Adam(self.parameters(),lr=lr, weight_decay=weight_decay)

        features= self.gc(torch.FloatTensor(X),torch.FloatTensor(adj))
        #----------------------------------------------------------------        
        if init=="kmeans":
            # print("Initializing cluster centers with kmeans, n_clusters known")
            self.n_clusters=n_clusters
            kmeans = KMeans(self.n_clusters, n_init=20)
            if init_spa:
                #------Kmeans use exp and spatial
                y_pred = kmeans.fit_predict(features.detach().numpy())
            else:
                #------Kmeans only use exp info, no spatial
                y_pred = kmeans.fit_predict(X)  #Here we use X as numpy
        elif init=="leiden":
            # print("Initializing cluster centers with leiden, resolution = ", res)
            if init_spa:
                adata=sc.AnnData(features.detach().numpy())
            else:
                adata=sc.AnnData(X)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
            sc.tl.leiden(adata,resolution=res)
            y_pred=adata.obs['leiden'].astype(int).to_numpy()
            self.n_clusters=len(np.unique(y_pred))
        #----------------------------------------------------------------
        y_pred_last = y_pred
        self.mu = Parameter(torch.Tensor(self.n_clusters, self.nhid))
        X=torch.FloatTensor(X)
        adj=torch.FloatTensor(adj)
        self.trajectory.append(y_pred)
        features=pd.DataFrame(features.detach().numpy(),index=np.arange(0,features.shape[0]))
        Group=pd.Series(y_pred,index=np.arange(0,features.shape[0]),name="Group")
        Mergefeature=pd.concat([features,Group],axis=1)
        cluster_centers=np.asarray(Mergefeature.groupby("Group").mean())
        
        self.mu.data.copy_(torch.Tensor(cluster_centers))
        self.train()
        for epoch in range(max_epochs):
            if epoch%update_interval == 0:
                _, q = self.forward(X,adj)
                p = self.target_distribution(q).data
            # if epoch%10==0:
            #     print("Epoch ", epoch) 
            optimizer.zero_grad()
            z,q = self(X, adj)
            loss = self.loss_function(p, q)
            loss.backward()
            optimizer.step()
            if epoch%trajectory_interval == 0:
                self.trajectory.append(torch.argmax(q, dim=1).data.cpu().numpy())

            #Check stop criterion
            y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
            delta_label = np.sum(y_pred != y_pred_last).astype(np.float32) / X.shape[0]
            y_pred_last = y_pred
            if epoch>0 and (epoch-1)%update_interval == 0 and delta_label < tol:
                # print('delta_label ', delta_label, '< tol ', tol)
                # print("Reach tolerance threshold. Stopping training.")
                # print("Total epoch:", epoch)
                break


    def fit_with_init(self, X,adj, init_y, lr=0.001, max_epochs=5000, update_interval=1, weight_decay=5e-4,opt="sgd"):
        # print("Initializing cluster centers with kmeans.")
        if opt=="sgd":
            optimizer = optim.SGD(self.parameters(), lr=lr, momentum=0.9)
        elif opt=="admin":
            optimizer = optim.Adam(self.parameters(),lr=lr, weight_decay=weight_decay)
        X=torch.FloatTensor(X)
        adj=torch.FloatTensor(adj)
        features, _ = self.forward(X,adj)
        features=pd.DataFrame(features.detach().numpy(),index=np.arange(0,features.shape[0]))
        Group=pd.Series(init_y,index=np.arange(0,features.shape[0]),name="Group")
        Mergefeature=pd.concat([features,Group],axis=1)
        cluster_centers=np.asarray(Mergefeature.groupby("Group").mean())
        self.mu.data.copy_(torch.Tensor(cluster_centers))
        self.train()
        for epoch in range(max_epochs):
            if epoch%update_interval == 0:
                _, q = self.forward(torch.FloatTensor(X),torch.FloatTensor(adj))
                p = self.target_distribution(q).data
            X=torch.FloatTensor(X)
            adj=torch.FloatTensor(adj)
            optimizer.zero_grad()
            z,q = self(X, adj)
            loss = self.loss_function(p, q)
            loss.backward()
            optimizer.step()

    def predict(self, X, adj):
        z,q = self(torch.FloatTensor(X),torch.FloatTensor(adj))
        return z, q




class GC_DEC(nn.Module):
    def __init__(self, nfeat, nhid1,nhid2, n_clusters=None, dropout=0.5,alpha=0.2):
        super(GC_DEC, self).__init__()

        self.gc1 = GraphConvolution(nfeat, nhid1)
        self.gc2 = GraphConvolution(nhid1, nhid2)
        self.dropout = dropout
        self.mu = Parameter(torch.Tensor(n_clusters, nhid2))
        self.n_clusters=n_clusters
        self.alpha=alpha

    def forward(self, x, adj):
        x=self.gc1(x, adj)
        x = F.relu(x)
        x = F.dropout(x, self.dropout, training=True)
        x = self.gc2(x, adj)
        q = 1.0 / ((1.0 + torch.sum((x.unsqueeze(1) - self.mu)**2, dim=2) / self.alpha) + 1e-6)
        q = q**(self.alpha+1.0)/2.0
        q = q / torch.sum(q, dim=1, keepdim=True)
        return x, q

    def loss_function(self, p, q):
        def kld(target, pred):
            return torch.mean(torch.sum(target*torch.log(target/(pred+1e-6)), dim=1))
        loss = kld(p, q)
        return loss

    def target_distribution(self, q):
        #weight = q ** 2 / q.sum(0)
        #return torch.transpose((torch.transpose(weight,0,1) / weight.sum(1)),0,1)e
        p = q**2 / torch.sum(q, dim=0)
        p = p / torch.sum(p, dim=1, keepdim=True)
        return p

    def fit(self, X,adj, lr=0.001, max_epochs=10, update_interval=5, weight_decay=5e-4,opt="sgd",init="leiden",n_neighbors=10,res=0.4):
        self.trajectory=[]
        # print("Initializing cluster centers with kmeans.")
        if opt=="sgd":
            optimizer = optim.SGD(self.parameters(), lr=lr, momentum=0.9)
        elif opt=="admin":
            optimizer = optim.Adam(self.parameters(),lr=lr, weight_decay=weight_decay)

        features, _ = self.forward(torch.FloatTensor(X),torch.FloatTensor(adj))
        #----------------------------------------------------------------
        
        if init=="kmeans":
            #Kmeans only use exp info, no spatial
            #kmeans = KMeans(self.n_clusters, n_init=20)
            #y_pred = kmeans.fit_predict(X)  #Here we use X as numpy
            #Kmeans use exp and spatial
            kmeans = KMeans(self.n_clusters, n_init=20)
            y_pred = kmeans.fit_predict(features.detach().numpy())
        elif init=="leiden":
            adata=sc.AnnData(features.detach().numpy())
            sc.pp.neighbors(adata, n_neighbors=n_neighbors)
            sc.tl.leiden(adata,resolution=res)
            y_pred=adata.obs['leiden'].astype(int).to_numpy()
        #----------------------------------------------------------------
        X=torch.FloatTensor(X)
        adj=torch.FloatTensor(adj)
        self.trajectory.append(y_pred)
        features=pd.DataFrame(features.detach().numpy(),index=np.arange(0,features.shape[0]))
        Group=pd.Series(y_pred,index=np.arange(0,features.shape[0]),name="Group")
        Mergefeature=pd.concat([features,Group],axis=1)
        cluster_centers=np.asarray(Mergefeature.groupby("Group").mean())
        
        self.mu.data.copy_(torch.Tensor(cluster_centers))
        self.train()
        for epoch in range(max_epochs):
            if epoch%update_interval == 0:
                _, q = self.forward(X,adj)
                p = self.target_distribution(q).data
            if epoch%100==0:
                print("Epoch ", epoch) 
            optimizer.zero_grad()
            z,q = self(X, adj)
            loss = self.loss_function(p, q)
            loss.backward()
            optimizer.step()
            self.trajectory.append(torch.argmax(q, dim=1).data.cpu().numpy())

    def fit_with_init(self, X,adj, init_y, lr=0.001, max_epochs=10, update_interval=1, weight_decay=5e-4,opt="sgd"):
        # print("Initializing cluster centers with kmeans.")
        if opt=="sgd":
            optimizer = optim.SGD(self.parameters(), lr=lr, momentum=0.9)
        elif opt=="admin":
            optimizer = optim.Adam(self.parameters(),lr=lr, weight_decay=weight_decay)
        X=torch.FloatTensor(X)
        adj=torch.FloatTensor(adj)
        features, _ = self.forward(X,adj)
        features=pd.DataFrame(features.detach().numpy(),index=np.arange(0,features.shape[0]))
        Group=pd.Series(init_y,index=np.arange(0,features.shape[0]),name="Group")
        Mergefeature=pd.concat([features,Group],axis=1)
        cluster_centers=np.asarray(Mergefeature.groupby("Group").mean())
        self.mu.data.copy_(torch.Tensor(cluster_centers))
        self.train()
        for epoch in range(max_epochs):
            if epoch%update_interval == 0:
                _, q = self.forward(torch.FloatTensor(X),torch.FloatTensor(adj))
                p = self.target_distribution(q).data
            X=torch.FloatTensor(X)
            adj=torch.FloatTensor(adj)
            optimizer.zero_grad()
            z,q = self(X, adj)
            loss = self.loss_function(p, q)
            loss.backward()
            optimizer.step()

    def predict(self, X, adj):
        z,q = self(torch.FloatTensor(X),torch.FloatTensor(adj))
        return z, q

class SpaGCN(object):
    def __init__(self):
        super(SpaGCN, self).__init__()
        self.l=None

    def set_l(self, l):
        self.l=l

    def train(self,adata,adj, embed,
            num_pcs=50, 
            lr=0.005,
            max_epochs=2000,
            weight_decay=0,
            opt="admin",
            init_spa=True,
            init="leiden", #leiden or kmeans
            n_neighbors=10, #for leiden
            n_clusters=None, #for kmeans
            res=0.4, #for leiden
            tol=1e-3):
        self.num_pcs=num_pcs
        self.res=res
        self.lr=lr
        self.max_epochs=max_epochs
        self.weight_decay=weight_decay
        self.opt=opt
        self.init_spa=init_spa
        self.init=init
        self.n_neighbors=n_neighbors
        self.n_clusters=n_clusters
        self.res=res
        self.tol=tol
        assert adata.shape[0]==adj.shape[0]==adj.shape[1]
        # pca = PCA(n_components=self.num_pcs)
        # if issparse(adata.X):
        #     pca.fit(adata.X.A)
        #     embed=pca.transform(adata.X.A)
        # else:
        #     pca.fit(adata.X)
        #     embed=pca.transform(adata.X)
        ###------------------------------------------###
        if self.l is None:
            raise ValueError('l should not be set before fitting the model!')
        adj_exp=np.exp(-1*(adj**2)/(2*(self.l**2)))
        #----------Train model----------
        self.model=simple_GC_DEC(embed.shape[1],embed.shape[1])
        self.model.fit(embed,adj_exp,lr=self.lr,max_epochs=self.max_epochs,weight_decay=self.weight_decay,opt=self.opt,init_spa=self.init_spa,init=self.init,n_neighbors=self.n_neighbors,n_clusters=self.n_clusters,res=self.res, tol=self.tol)
        self.embed=embed
        self.adj_exp=adj_exp

    def predict(self):
        z,q=self.model.predict(self.embed,self.adj_exp)
        y_pred = torch.argmax(q, dim=1).data.cpu().numpy()
        # Max probability plot
        prob=q.detach().numpy()
        return y_pred, prob