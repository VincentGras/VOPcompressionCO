# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 10:27:00 2022

@author: VG244206
"""

import numpy as np
import nlopt
from os.path import join as pjoin
import scipy.io as sio
import multiprocessing as mp
import time

def isquaremat(Q, Nc):
    return  Q.ndim == 2 and Q.shape[0]==Nc and Q.shape[1]==Nc

def computeVOP(pool, criterion, lQ, Qmargin, c0 = None, Qvop = None):
    
    # compute VOP using convex optimization (CO)
    # input :
    #    pool : multiprocessing object
    #    lQ :    list (length N) of NcxNc Hermitian matrices :   the 10g-SAR matrices
    #    Qmargin : NcxNc Hermitian matrix  : compression margin
    #    c0 :   list of integers  : initial classification (empty by default)
    #    lQvop
    # output :
    #    c   :   list (length N) of integers : output classification
    
    # Number of SAR matrices
    N = len(lQ);
    
    if N==0:
        c = []
        return
    
    # Matrix Dimension
    Nc = lQ[0].shape[0]
    
    for Q in lQ:
        assert isquaremat(Q, Nc), "lQ : expecting list of square {Nc}x{Nc} matrices".format(Nc=Nc);
               
    assert isquaremat(Qmargin, Nc), "Qmargin : expecting a square {Nc}x{Nc} matrices".format(Nc=Nc);          
        
    # Classification
    
    nyc = np.uint8(2)
    dominated = np.uint8(0)
    vop = np.uint8(1)
    
    if (not (c is None)):
        assert c.size==N, 'Bad initial classification'
    else:
        c = np.ndarray((N))
        c[:] = nyc
           
    if (not (Qvop is None)):
        for Q in Qvop:
            assert isquaremat(Q, Nc), "Qvop : expecting list of square {Nc}x{Nc} matrices".format(Nc=Nc);   
            
    remain = N
    
    while (remain>0) :
        
        print('%d new VOPs / %d Q-matrices / %d to test still' % (np.sum(c==vop), N, remain));
        
        # List SAR matrices that are not yet classified 
        J = list(np.nonzero(c==nyc)[0]);
        
        if (Qvop is None):
            R = pool.starmap(criterion, [(lQ[j],[Qmargin]) for j in J]);
        else:
            R = pool.starmap(criterion, [(lQ[j],Qvop) for j in J]);
        
        rmax = 0;
        jmax = -1;
        
        for j,r in zip(J,R):
            
            if (r<=1):
                
                c[j]=dominated;
                remain = remain - 1;
                
            elif r>rmax:
                
                jmax=j;
                rmax=r;
                
        # if rmax>0, mark as new vop the SAR matrix that realizes the max 
        
        if (rmax>0):
            
            c[jmax] = vop; 
            Qvop.append (lQ[jmax]) + Qmargin;
            remain = remain - 1;
       
        
    print('%d new VOPs / %d Q-matrices' % (np.sum(c==vop), N));        
    
    return c, Qvop   

def SARC(Q, lVOP):

    # Test if Q <= Qvop (combined domination)
    # input :
    #    Q10g : NcxNc Hermitian matric 
    #    Qvop : list of Nvop NcxNc Hermitian matrices
    # Output :
    #    R       N real>0      max ((V^H Q V)/(V^H Qvop V))
    #    V       NcxN complex    argmax ((V^H Q V)/(V^H Qvop V))

    
    Nvop = len(lVOP);            
    Nc = Q.shape[0];
    
    # Convert to real symetric matrices
    CHtoRS = lambda Q: np.vstack((np.hstack((np.real(Q), -np.imag(Q))), np.hstack((np.imag(Q), np.real(Q)))));
    
    QvopR = [CHtoRS(Qvop) for Qvop in lVOP];
    QR = CHtoRS(Q);
    
    SAR = lambda Qvop,V: V.transpose() @ Qvop @ V;
        
    if  Nvop == 1:
            
        # in this case do an eignevalue analysis of QR-QvopR
        # and look for the maximum eigenvalue
        # If that eigenvalue is negative then QR<=QvopR and R[i] will be < 1
         
        H = QR-QvopR[0];
        eigW,eigV = np.linalg.eigh(H);
        VR = eigV[:,2*Nc-1]; # vector corresponding to the max eigenvalue
        
    else:
            
        # In this case run a convex optimization
        
        V0R = np.zeros((2*Nc));
        V0R[0] = 1;
        eigW,eigV = np.linalg.eigh(QR);
        V0R = eigV[:,2*Nc-1];
        C = max([SAR(Qvop,V0R) for Qvop in lVOP]);
        V0R *= 0.5 / (np.sqrt(max(0,C))+1e-9);
        
        # prepare nlopt     
        opt = nlopt.opt(nlopt.LD_MMA , 2*Nc)
        opt.set_min_objective(lambda x,grad: objfun(QR, x, grad))
        opt.add_inequality_mconstraint(lambda res,x,grad: constrfun(QvopR, res, x, grad), 1e-8*np.ones(Nvop))
        opt.set_xtol_rel(1e-6)
        opt.set_maxeval(200)
        
        # find minimizer
        VR = opt.optimize(V0R)
        
    # Store solution of the problem  max ((V^H Q V)/(V^H Qvop V))
    R = SAR(QR,VR) / SAR(QvopR,VR) - 1;
    V = VR[0:Nc] + 1j * VR[Nc:2*Nc];

    return (R,V)

# Objective function

def objfun(Q, V, G) :
    P = -V.transpose() @ Q
    C = P @ V;
    G[:] = 2 * P  
    
    return C;


# Constraints

def constrfun(lQ,C,V,G) :
    
    Vt = V.transpose();
    
    for i,Q in zip(range(0,len(lQ)),lQ) :
        
        P = Vt @ lQ[i];
        C[i] = P @ V - 1;
        G[i,:] = 2*P;


# main : test code in a simple example

if __name__ == '__main__': 
    
    
    data_dir_QtoTest = '.'
    mat_fname = pjoin(data_dir_QtoTest, 'test_Q.mat')
    mat_contents = sio.loadmat(mat_fname)
    Q = mat_contents['Q10g']
    Q = Q.transpose((2, 0, 1))
    Q = list(Q)
    
    Nc = Q[0].shape[0]
   
    Qm = (1/250)*np.eye(Nc,dtype='complex64')
    
    
    # time 
    t = time.time()
    
    # open mp pool
    pool = mp.Pool(mp.cpu_count())
    print("nb of procs = ", mp.cpu_count())
        
    vop, Qvop = computeVOP(pool,SARC,Q,Qm)
    
    # close pool
    pool.close()
    
      # elapsed time
    elapsed = time.time() - t
    print('elapsed time VOP compression) ' + str(elapsed) + ' s')
    
    
   
    
    
    
