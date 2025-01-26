"""
sinkhorn_log - stabilized sinkhorn over log domain with acceleration
[u,v,err] = sinkhorn_log(mu,nu,c,epsilon,options);
mu and nu are marginals.
c is cost
epsilon is regularization
coupling is 
gamma = exp( (-c+u*ones(1,N(2))+ones(N(1),1)*v')/epsilon );

options.niter is the number of iterations.
options.tau is an avering step. 
  - tau=0 is usual sinkhorn
  - tau<0 produces extrapolation and can usually accelerate.
options.rho controls the amount of mass variation. Large value of rho
impose strong constraint on mass conservation. rho=Inf (default)
corresponds to the usual OT (balanced setting). 

Copyright (c) 2016 Gabriel Peyre

@author: linxin,caiqy
"""
import numpy as np
from multiprocessing import Pool
from scipy.linalg.blas import sgemm
def sinkhorn_log(mu, nu, c, epsilon=1e-6, niter=50):
    lamb = 1e4/(1e4+1) # 1e4 is the default value in the original code
    N = np.asarray([mu.size, nu.size])
    H = np.ones((N[0],N[1]))
    u = np.zeros(N[0])
    v = np.zeros(N[1])
    for i in range(niter):
        M = (-c + fast_dot(np.diag(u) , H) + fast_dot(H , np.diag(v)))/epsilon
        u = ave(u, lamb*epsilon*np.log(mu+1e-20) - lamb*epsilon*lse(M) + lamb*u, 1/2)
        v = ave(v, lamb*epsilon*np.log(nu+1e-20) - lamb*epsilon*lse(M.transpose()) + lamb*v, 1/2)
    gamma = np.exp(M)
    return gamma

def ave(u, u1, tau):
    return tau*u+(1-tau)*u1

def lse(A):
    return np.log(np.exp(A).sum(axis=1)+1e-100)

def fast_dot(a,b):
    return sgemm(1, a = a, b = b)

def _solve_ot(k):
    return np.array(sinkhorn_log(scy[k,:], scx[k,:], grid[:,:], 1e-6, 50))

def solve_OT(emptyGrids,scLocX,scLocY,thres=0.01,njob=8):
    """
    This function is used to solve the OT problem.
    The output is a matrix with the same shape as the spatial transcriptomics data, where each element is a probability of the corresponding cell belonging to a spoy.
    """
    # grid需要增加同个权重还是不同权重的区别
    global grid
    grid = emptyGrids
    scLocX_0 = scLocX*(scLocX>thres).astype(int)
    scLocY_0 = scLocY*(scLocY>thres).astype(int)
    zx = scLocX_0.sum(axis=1)
    zy = scLocY_0.sum(axis=1)
    global scx,scy
    scLocX_sel = scLocX[(zx!=0) & (zy!=0),]
    scLocY_sel = scLocY[(zx!=0) & (zy!=0),]
    scx = scLocX_sel * ((scLocX_sel > 0).astype(int))
    scy = scLocY_sel * ((scLocY_sel > 0).astype(int))
    print('Starting solve the OT problem.')
    # multiprocessing
    pool = Pool(njob) 
    gamma = pool.map(_solve_ot, range(scLocX_sel.shape[0]))
    pool.close()
    pool.join()
    gamma = np.array(gamma)
    return gamma