from os.path import dirname, join as pjoin
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as pypl
import time
from  computeVOP_CO import computeVOP_CO__mp
import multiprocessing as mp
    

if __name__ == '__main__': 
    
    
    # générer Q10gToTest
    
    
    data_dir_QtoTest = '.'
    mat_fname = pjoin(data_dir_QtoTest, 'Ella_Q.mat')
    mat_contents = sio.loadmat(mat_fname);
    Q10g = mat_contents['Q10g']
    
    Nc = Q10g.shape[0]
   
    Qm = (1/250)*np.eye(8,dtype='complex64');
    
    
    # time 
    t = time.time();
    
    # open mp pool
    pool = mp.Pool(mp.cpu_count())
    print("nb of procs = ", mp.cpu_count());
    # run VOp compression using mp
    Qvop = computeVOP_CO__mp(pool,Q10g[:,:,range(0,Q10g.shape[2])],Qm);
    # close pool
    pool.close();
    
      # elapsed time
    elapsed = time.time() - t
    print('elapsed time VOP compression) ' + str(elapsed) + ' s')
    
    
    # verif
    
    Q10gL = [];
    for i in range(0,Q10g.shape[2]):
        Q10gL.append(Q10g[:,:,i]);

    # time 
    t = time.time();
    
    RV = pool.starmap(testQmatrixDomination_CO, [(Q,Qvop) for Q in Q10gL])
    
    # Step 3: Don't forget to close
    pool.close()    

    # elapsed time verif
    elapsed = time.time() - t
    print('elapsed time VOP verif) ' + str(elapsed) + ' s')
    # viz
    
    R = np.zeros(len(RV))
    i=0;
    for rv in RV:
        R[i]=rv[0];
        i+=1


    pypl.plot(R)
    
    
    