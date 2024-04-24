'''
Filename    : LLE by LOBPCG.py
Description : implementation of locally linear embedding and LOBPCG
Date        : 2022.1.16
'''
import numpy as np
from sklearn import neighbors
from scipy.linalg import solve, eigh
from scipy.sparse import linalg
from scipy import sparse
from time import time

######################################
#                LOBPCG
######################################

def Orthonormalize(A):
    return np.linalg.qr(A)[0]


def LOBPCGsvd(A,x,maxiter = 1000,retFrobeniusHistory=False):
    '''
    Find the smallest x singular vectors and singular values of A

    Parameters
    ----------
    A : a scipy sparse matrix (Warning : numpy array is forbidden)

    x : target dimension or initial guess of singular vectors

    maxiter : maximum iteration time

    retFrobeniusHistory : whether or not return the frobenius norm of Ax

    Returns 
    -------
    x : the eigenvectors

    mu : corresponding eigenvalues

    history : returns the Frobenius history iff retFrobeniusHistory == True
    '''

    n = A.shape[0]
    if type(x) == int:
        m = x
        x = np.random.randn(n*m).reshape((n,m))
    else:
        m = x.shape[1]
        x = x.copy()
    x = Orthonormalize(x)
    
    x2 = np.linalg.eigh((x.T * A.T) @ (A * x))
    mu = x2[0]
    x2 = x2[1]
    x = x @ x2

    p = np.zeros((n,m))

    history = []
    for k in range(maxiter):
        Ax = A * x
        r = A.T * Ax - mu * x  # broadcast
        #print('mu =',mu)
        #print('re =',np.linalg.norm(r[:,0]))#,np.linalg.norm(r[:,1])
        r -= x @ (x.T @ r)
        #r = Tr  # preconditioner
        r = Orthonormalize(r)

        if k != 0:
            Ap = A * p
            if retFrobeniusHistory:
                history.append(np.linalg.norm(Ax))
            Ar = A * r
            if k == 1:
                gramA = np.zeros((3*m,3*m))
            gramA[:m,:m].flat[::m+1] = mu   # = np.diag(mu)
            gramA[:m,m:-m] = Ax.T @ Ar
            gramA[:m,-m:] = Ax.T @ Ap
            gramA[m:-m,m:-m] = Ar.T @ Ar
            gramA[m:-m,-m:] = Ar.T @ Ap
            gramA[-m:,-m:] = Ap.T @ Ap

            gramA[m:-m,:m]  = gramA[:m,m:-m].T
            gramA[-m:,:m]   = gramA[:m,-m:].T
            gramA[-m:,m:-m] = gramA[m:-m,-m:].T


            gramB = np.eye(3*m)
            gramB[:m,-m:] = x.T @ p
            gramB[m:-m,-m:] = r.T @ p

            gramB[-m:,:m]   = gramB[:m,-m:].T
            gramB[-m:,m:-m] = gramB[m:-m,-m:].T
        else: # k == 0:
            Ar = A * r

            gramA = np.zeros((2*m,2*m))
            gramA[:m,:m] = np.diag(mu)
            gramA[:m,m:] = Ax.T @ Ar
            gramA[m:,m:] = Ar.T @ Ar

            gramA[m:,:m] = gramA[:m,m:].T


            gramB = np.eye(2*m)
            gramB[:m,m:] = x.T @ r
            
            gramB[m:,:m] = gramB[:m,m:].T


        x2 = eigh(gramA,gramB,check_finite=False)
        #print(x2)
        iid = np.argsort(x2[0])[:m] 
        mu = x2[0][iid]
        x2 = x2[1][:,iid]

        if k != 0:
            p = p @ x2[-m:]
            p += r @ x2[m: -m]
        else: # k == 0:
            p = r @ x2[m:]
        
        x = x @ x2[:m] + p
        p = Orthonormalize(p)
        
    if retFrobeniusHistory:
        return x, mu, history
    return x, mu
        

######################################
#                 llE
######################################

def KNearest(X,k):
    nbors = neighbors.NearestNeighbors(n_neighbors=k+1,algorithm='ball_tree').fit(X)
    nbors = nbors.kneighbors(X,return_distance=False)[:,1:]
    return nbors

def ComputeWeight(X,Z,t,overwrite):
    C = X[Z] - X[t]
    G = C @ C.T
    trace = np.trace(G)
    reg = min(1e-3, 1e-3*trace)
    G.flat[::G.shape[0] + 1] += reg
    overwrite[:] = solve(G, np.ones(G.shape[0], dtype=G.dtype), sym_pos=True)
    overwrite *= 1. / np.sum(overwrite)

def GetWeights(X,k,Z):
    n = X.shape[0]
    w_data = np.zeros(n*k, dtype=X.dtype)
    for i, j in enumerate(range(0,n*k,k)):
        ComputeWeight(X,Z[i],i,w_data[j:j+k])
    return w_data

def Embedding(W,d,maxiter=1000,retFrobeniusHistory=False):
    #return np.linalg.svd(W.A)[2][-d-1:-1].T
    #eigs = linalg.lobpcg(W, np.random.randn(n*(d+1)).reshape((n,d+1)), largest=False)
    return LOBPCGsvd(W,d+1,maxiter=maxiter,retFrobeniusHistory=retFrobeniusHistory)

def LocallyLinearEmbedding(X,k,d,maxiter = -1,verbosityLevel = 1,
        retW = False,retFrobeniusHistory = False):
        
    '''
    Apply LLE on data X

    Parameters
    ----------
    X : data matrix

    k : neighborhood size

    d : embed dimension 

    maxiter : maximum iteration time

    verbosityLevel : whether or not show the progress

    retW : whether or not return the weights

    retFrobeniusHistory : whether or not return the history of LOBPCG

    Returns 
    -------
    y : the embedded data

    W : returns the weights iff retW == True

    history : returns the Frobenius history iff retFrobeniusHistory == True
    '''

    n = X.shape[0]

    if maxiter < 0:
        maxiter = int(1000 * max(1, np.log(n) / 6.9))


    if verbosityLevel:
        print('Searching for K-Neighbors ... ', end = '')
        clock = time()
    nbors = KNearest(X,k)

    if verbosityLevel:
        print(' %f s\nComputing weights ... ' % (time() - clock), end = '')
        clock = time()
        
    w_data = GetWeights(X,k,nbors)
    W = sparse.csr_matrix((w_data, nbors.flatten(), np.arange(0,(n+1)*k,k)) , shape=(n,n))
    W = sparse.eye(n) - W
    #return W
    
    if verbosityLevel:
        print(' %f s\nComputing embeddings ...' % (time() - clock), end = '')
        clock = time()
    embeds = Embedding(W,d,maxiter,retFrobeniusHistory=retFrobeniusHistory)
    embedvec = embeds[0] * np.sqrt(n)
    
    if verbosityLevel:
        print(' %f s' % (time() - clock))
    
    if retFrobeniusHistory:
        if retW:
            return embedvec, W, embeds[2]
        else:
            return embedvec, embeds[2]
       
    if retW:
        return embedvec, W
    return embedvec
