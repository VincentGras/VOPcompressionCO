from os.path import dirname, join as pjoin
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as pypl
import time
from  computeVOP_CO import computeVOP_CO__mp
import multiprocessing as mp
from  testQmatrixDomination_CO import testQmatrixDomination_CO


if __name__ == '__main__': 
    
    
    # générer Q10gToTest
    
    
    data_dir_QtoTest = '.'
    mat_fname = pjoin(data_dir_QtoTest, 'Ella_Q.mat')
    mat_contents = sio.loadmat(mat_fname);
    Q = mat_contents['Q10g']
    Q = Q.transpose((2, 0, 1))
    Q = list(Q);
    
    Nc = Q[0].shape[0]
   
    Qm = (1/250)*np.eye(Nc,dtype='complex64');
    
    
    # time 
    t = time.time();
    
    # open mp pool
    pool = mp.Pool(mp.cpu_count())
    print("nb of procs = ", mp.cpu_count());
    # run VOp compression using mp
    
    
    r = testQmatrixDomination_CO(Q[0],Qm)
    
    vop = computeVOP_CO__mp(pool,Q,Qm);
    
    # close pool
    pool.close();
    
      # elapsed time
    elapsed = time.time() - t
    print('elapsed time VOP compression) ' + str(elapsed) + ' s')
    
    
   
    
    
    