'''
Filename    : Schur Decomposition.py
Description : Computing the Schur Decomposition by Francis implicit double shift QR algorithm
Author      : Forever Haibara
Date        : 2021.12.31
'''

import numpy as np
from numba import jit
from math import sqrt
from sys import float_info
from tqdm import tqdm

@jit(nopython = True)
def squarenorm(x):
    return np.dot(x,x) # np.dot(x.conj(),x).real

@jit(nopython = True)
def cnorm(x):
    return x * x # x.real * x.real + x.imag * x.imag
    
#@jit(nopython = True,cache=True)
@jit(nopython=True)
def HouseholderVector(A):
    '''
    Compute real Householder vector of norm 2. 
    
    If no Householder is needed, the second argument returns one.
    '''
    squarelength = squarenorm(A[1:]) 
    if squarelength == 0:
        return A, 1 # A, A, 1
    normal     = A.copy()
    rotation   = 1 if A[0]>=0 else -1
    normal[0] = A[0] + sqrt(squarelength + cnorm(A[0]))*rotation
    normal *= sqrt(2. / (squarelength + cnorm(normal[0])))
    return normal, 0

@jit(nopython=True)
def leftHouse(A,normal):
    '''
    left Householder transformation where ||``normal``|| = sqrt 2
    '''
    tmp = normal @ A
    n = tmp.shape[0]
    for i in range(normal.shape[0]):
        for j in range(n):
            A[i,j] -= normal[i]*tmp[j]

@jit(nopython=True)
def rightHouse(A,normal):
    '''
    right Householder transformation where ||``normal``|| = sqrt 2
    '''
    tmp = A @ normal
    n = normal.shape[0]
    for i in range(tmp.shape[0]):
        for j in range(n):
            A[i,j] -= tmp[i]*normal[j]

#@jit(nopython = True,cache=True)
@jit(nopython=True)
def Hessenberg(A):
    ''' 
    Convert a matrix to upper-Hessenberg
    '''
    n = A.shape[0] 
    Q = np.eye(n,dtype=A.dtype)
    for i in range(n-1):
        normal, invalid = HouseholderVector(A[i+1:,i])
        if invalid:
            continue
        # A ->  [I- (2nn*)] A   =   A - 2n(n* A)
        leftHouse(A[i+1:,i:],normal)
        A[i+2:,i] = 0
        rightHouse(A[:,i+1:],normal)
        rightHouse(Q[:,i+1:],normal)
    return (Q,A)
    
@jit(nopython=True)
def ImplicitDoubleShift(A,Q,U,V):
    '''
    Perform a double shift on matrix A and update A,Q,U,V
    '''
    s = A[-1,-1] + A[-2,-2]
    t = A[-1,-1]*A[-2,-2] - A[-1,-2]*A[-2,-1]
    normal, invalid = HouseholderVector(np.array([
                A[0,0]*(A[0,0]-s) + A[0,1]*A[1,0] + t,
                A[1,0]*(A[0,0]+A[1,1]-s),
                A[1,0]*A[2,1]
                        ], dtype=A.dtype))
    if invalid:
        return
    leftHouse(A[:3,:],normal)
    rightHouse(A[:4,:3],normal)
    leftHouse(U[:3,:],normal)
    rightHouse(V[:,:3],normal)
    rightHouse(Q[:,:3],normal)

@jit(nopython=True)
def BuldgeChasing(A,Q,U,V):
    '''
    Convert the matrix back to Hessenberg and update A,Q,U,V
    '''
    n = A.shape[0]
    for i in range(1,n-1): 
        j = i + 3
        normal, invalid = HouseholderVector(A[i:j,i-1])
        if invalid:
            continue

        leftHouse(A[i:j,:],normal)
        A[i+1:j,i-1] = 0
        rightHouse(A[:j+1,i:j],normal)
        leftHouse(U[i:j,:],normal)
        rightHouse(V[:,i:j],normal)
        rightHouse(Q[:,i:j],normal)
        

@jit(nopython = True)
def Deflation(A,start,end,eps):
    '''
    Determine the eigenvalues to deflat and an irreduced Hessenberg matrix to operate on.
    '''
    for i in range(start,end):
        if A[i,i-1] != 0. and abs(A[i,i-1]) < eps * (abs(A[i,i])+abs(A[i-1,i-1])):
            A[i,i-1] = 0.
        
    for i in range(end-1,0,-1):
        if A[i,i-1] == 0.:
            if end - i <= 2:
                end = i
            else :
                start = i
                break
    else :
        start = 0
    return start , end

#@jit
def SchurDecomposition(A,maxiter=-1,copy=1,verbose=-1,printepoch=0):
    '''
    Compute Schur Decomposition for a real matrix. Return Q,T such that Q*AQ = T.

    ``maxiter`` the maximum iteration steps permitted, defaults to 3n + 60

    ``copy`` overwrites on A if ``copy`` == 0, defaults to 1

    ``verbose`` verbose the progress if True, defaults to True when n >= 500. 

    ``printepoch`` show how many iterations it takes if True, defaults to False
    '''
    n , lastend = A.shape
    start , end = 0 , lastend

    if copy: A = A.copy()
    if verbose < 0: verbose = 1 if n >= 500 else 0
    if maxiter < 0: maxiter = 3*n + 60
    eps = float_info.epsilon
    Q , A = Hessenberg(A)
    i = 0

    if verbose:
        pbar =  tqdm(total = n)
        pbar.set_description('Deflation')

    for i in range(maxiter):
        if end <= 2:
            break

        ImplicitDoubleShift(A[start:end,start:end],Q[:,start:end],
                            A[start:end,end:],A[:start,start:end])
        BuldgeChasing(A[start:end,start:end],Q[:,start:end],
                            A[start:end,end:],A[:start,start:end])

        # deflation
        start , end = Deflation(A,start,end,eps)

        if verbose and i % 10 == 0:
            pbar.set_postfix(steps = i)
            pbar.update(lastend - end)
            lastend = end
            
    else :
        print('Warning: Convergence Failure')
        
    if verbose:
        pbar.update(lastend)
        pbar.close()
    
    if printepoch:
        print('Steps =',i)

    return (Q , A)

# Sample
if __name__ == '__main__':
    n = 500
    np.random.seed(0)
    A = np.random.randn(n*n).reshape((n,n))
    Q , T = SchurDecomposition(A,verbose=1,printepoch=1)
    print('Forward Error =', np.linalg.norm(Q.conj().T @ A @ Q - T))
    print('Orthogonality Loss =',np.linalg.norm(Q.T @ Q - np.eye(n)))