from os.path import dirname, join as pjoin
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import time
from  VOPcompressionCO import compute_rqstar, compute_vop
import multiprocessing as mp


if __name__ == '__main__':


    # générer Q10gToTest


    data_dir_QtoTest = '.'
    mat_fname = pjoin(data_dir_QtoTest, 'testQ.mat')
    mat_contents = sio.loadmat(mat_fname)
    Q = mat_contents['Q10g']
    lQ = list(Q.transpose((2, 0, 1)))

    Nc = lQ[0].shape[0]
    print(Nc)

    Qm = (1 / 250) * np.eye(Nc, dtype='complex128');

    pool = None
    # open mp pool
    # pool = mp.Pool(2)
    # print("nb of procs = ", mp.cpu_count());
    lVOP, classif, rqstar = compute_vop(lQ, compression_margin=Qm, pool=pool)
    # pool.close()
    r = [compute_rqstar(_Q, lVOP)[0] for _Q in lQ]
    plt.figure()
    plt.plot(np.array(r))

    # lVOP = []
    # updateVOP(lVOP, lQ, R, Qmargin=Qm, maxR=1.1)
    # print(len(lVOP))
