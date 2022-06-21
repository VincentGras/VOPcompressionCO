import numpy as np
from  computeVOP_CO import computeVOP_CO__mp


def computeVOP_iCO(Q10g, Qmargin, niter, r) :
 
 
    Nc = Q10g.shape[0];
    
    # nb of SAR matrices
    
    if Q10g.ndim<3:
        Q10g = Q10g.reshape((Nc,Nc,1));
    N = Q10g.shape[2];
    
    # intialize Qvop to []
    Qvop = np.ndarray((Nc,Nc,0));
    
    # initialize margin
    Qm = (r ** (1-niter)) * Qmargin;
    
    # classif
    c = np.zeros((N));
  
    for i in range(1,niter) :
        
        
        Qvop, c = computeVOP_CO__mp(Q10g, Qmargin, Q10g[:,:,c>0] + Qm);
        print ('Nb of VOP after iteration %d/%d = %d' % i, niter, np.sum(c))
        Qm *= r;
        
