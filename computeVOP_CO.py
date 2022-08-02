# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 10:27:00 2022

@author: VG244206
"""

import numpy as np
from  testQmatrixDomination_CO import testQmatrixDomination_CO

def computeVOP_CO__mp(pool, lQ, Qmargin, c0=np.ndarray((0),dtype='uint8')):
    
    # compute VOP using convex optimization (CO)
    # input :
    #    pool : multiprocessing object
    #    lQ :    list (length N) of NcxNc Hermitian matrices :   the 10g-SAR matrices
    #    Qmargin : NcxNc Hermitian matrix  : compression margin
    #    c0 :   list of integers  : initial classification (empty by default)
    # output :
    #    c   :   list (length N) of integers : output classification
    
    # Number of SAR matrices
    N = len(lQ);
    
    if N==0:
        c = []
        return
    
    # Matrix Dimension
    Nc = lQ[0].shape[0];
    
    for Q in lQ:
        assert Q.ndim == 2 and Q.shape[0]==Nc and Q.shape[1]==Nc, 'lQ : expecting list of square matrices'    
               
    # Classification
    c = np.zeros(N, dtype='uint8');
    nyc = np.uint8(2);
    dominated = np.uint8(0);
    vop = np.uint8(1);
    
    if (c0.size>0):
        assert c0.size==N, 'Bad c0 (does not match Q10g 3rd dim)'
        c[:] = c0[:];
    else:
        c[:] = nyc;
    
           
    # initialize Qvop
    Qvop = [];
    for Q in lQ:
        Qvop.append(Q+Qmargin);

    
    remain = N;
    
    while (remain>0) :
        
        print('%d new VOPs / %d Q-matrices / %d to test still' % (np.sum(c==vop), N, remain));
        
        # List SAR matrices that are not yet classified 
        J = list(np.nonzero(c==nyc)[0]);
        
        if (len(Qvop)==0):
            R = pool.starmap(testQmatrixDomination_CO, [(lQ[j],[Qmargin]) for j in J]);
        else:
            R = pool.starmap(testQmatrixDomination_CO, [(lQ[j],Qvop) for j in J]);
        
        rmax=0;
        jmax=-1;
        
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
    
    return c   

