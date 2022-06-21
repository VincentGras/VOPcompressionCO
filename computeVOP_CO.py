# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 10:27:00 2022

@author: VG244206
"""

import numpy as np
from  testQmatrixDomination_CO import testQmatrixDomination_CO

def computeVOP_CO__mp(pool, Q10g, Qmargin, Qvop=np.ndarray((0))):
    
    # compute VOP using convex optimization (CO)
    # input :
    #    pool : multiprocessing object
    #    Q10g : NcxNxN ndarray   the 10g-SAR matrices
    #    Qmargin : NxxNc ndarray  compression margin
    #    Qvop : NcxNcxNvopInit initial VOPs (empty by default)
    # output :
    #   Qvop : NcxNcxNvop ndarray  the VOPs
    #   c    : SAR matrix classification (c(i)=1 if Q(:,:,i) is a VOP, 0 otherwise)
    #   Note that each added VOP is of the form Q(:,:,i)+Qmargin
    
    import numpy as np
    
    # Matric Dimension
    Nc = Q10g.shape[0];
    assert Q10g.shape[1]==Nc, 'Q10g : expecting square matrices'    
    
    # Number of SA matrices
    N = 1;
    
    if (Q10g.ndim>2) :
        N = Q10g.shape[2];
    else :
        Q10g = Q10g.reshape((Nc, Nc, 1));
        
    # Classification
    c = np.zeros(N, 'int8');
    nyc = 2;
    dominated = 0;
    vop = 1;
    c[:] = nyc;
    
    # R values that are computed on each SAR matrix for VOP decision
    R = np.zeros(N);
    R[:] = np.Inf;
    
    if (Qvop.size>0):
        
        assert Qvop.shape[0] == Nc, 'Bad Qvop (does not match Q10g in shape)'
        assert Qvop.shape[1] == Nc, 'Bad Qvop (does not match Q10g in shape)'
        assert Qvop.ndim <= 3, 'Bad Qvop shape (must be 2 or 3)'        
    
        if (Qvop.ndim==2):
            Qvop = Qvop.reshape((Nc,Nc,1));
    
       
    remain = N;
    
    while (remain>0) :
        
        print('%d new VOPs / %d Q-matrices / %d to test still' % (np.sum(c==vop), N, remain));
        
        # List Q10g matrices that are not yet classified 
        J = np.nonzero(c==nyc);
        Qs = [];
        if (Q10g.ndim>2) :
            for j in J[0]:
                Qs.append(Q10g[:,:,j]);
        else:
            Qs.append(Q10g);
        
        # run convex optimization
        if (Qvop.size <= 0):
            res_mp = pool.starmap(testQmatrixDomination_CO, [(Q,Qmargin) for Q in Qs]);
        else :
            res_mp = pool.starmap(testQmatrixDomination_CO, [(Q,Qvop) for Q in Qs]);
       
        # Classifiy
        rmax = 0;
        k = N+1;
        
        for (res, j) in zip(res_mp, J[0]):
            R[j] = res[0];
            
            if (R[j]>rmax):
                rmax=R[j];
                k=j;
            
            if (R[j] <= 1):
                c[j] = dominated;
                remain = remain - 1;
           
        if (rmax > 1) :
            c[k] = vop;
            remain = remain - 1;
        
        if Qvop.size == 0 :
            Qvop = Q10g[:,:,k] + Qmargin;
            Qvop = Qvop.reshape((Nc,Nc,1));
        else:
            newQvop = Q10g[:,:,k]+Qmargin;
            Qvop = np.concatenate( (newQvop.reshape((Nc,Nc,1)), Qvop), 2);
        
    
    print('%d new VOPs / %d Q-matrices' % (np.sum(c==vop), N));        
    
    return Qvop,c;   

