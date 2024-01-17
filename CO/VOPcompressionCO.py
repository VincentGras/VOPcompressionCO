# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 14:15:56 2023

@author: VG244206
"""

import numpy as np
import scipy.optimize
import numbers


def compute_vop(
    lQ,
    lVOP=[],
    real=None,
    classif=None,
    rqstar=None,
    spectral_norm=None,
    compression_margin=0.0,
    rqstar_margin=0.0,
    max_vop_nb=np.Inf,
    pool=None,
):
    # Number of SAR matrices
    N = len(lQ)
    if N == 0:
        return lVOP, np.zeros((0,), dtype=int), None
    # real symetric matrices ?
    if real is None:
        real = _is_real(lQ)
    elif real:
        assert _is_real(lQ), "expecting real symmetric SAR matrices"
    if real and len(lVOP) > 0:
        assert _is_real(lVOP), "expecting real symmetric VOP matrices"
    # Matrix dimension
    Nc = lQ[0].shape[0]
    # check consistency
    for (i, Q) in zip(range(N), lQ):
        assert (
            Q.ndim == 2 and Q.shape[0] == Nc and Q.shape[1] == Nc
        ), "lQ : expecting list of square matrices"
    # spectral norms
    if spectral_norm is None:
        spectral_norm = np.ndarray((N), dtype=float)
        for (i, Q) in zip(range(N), lQ):
            spectral_norm[i] = np.linalg.norm(Q, 2)
    else:
        assert isinstance(spectral_norm, np.ndarray), "expecting numpy array"
        assert spectral_norm.ndim == 1, "expecting 1d numpy array"
        assert spectral_norm.size == N, "expecting (%d,) numpy rray" % N
    max_spectral_norm = np.max(spectral_norm)
    # Compression margin
    if isinstance(compression_margin, np.ndarray):
        assert compression_margin.ndim == 2, "expecting matrix"
        assert (
            compression_margin.shape[0] == Nc
        ), "expecting a % d x % d matrix" % (Nc, Nc)
        assert (
            compression_margin.shape[1] == Nc
        ), "expecting a % d x % d matrix" % (Nc, Nc)
    else:
        assert isinstance(
            compression_margin, numbers.Number
        ), "expecting a scalar or a %d  x%d matrix" % (Nc, Nc)
        compression_margin = compression_margin * np.eye(Nc)
    # rqstar array
    if rqstar is None:
        rqstar = (np.zeros((N,), dtype=float),
                  np.zeros((N, Nc), dtype=lQ[0].dtype))
        # first element encodes the rq* value
        # 2nd to last element encodes the argmax (N-vector)
    # classif
    if classif is None:
        classif = np.zeros((N), dtype=int) - 1
    # Initialize indices of the remaining matrices
    remain = np.where(classif < 0)[0]
    # local compye_rqstar with only to positional args: Q and initialization

    # _compute_rqstar(Q, lVOP, real, scaling, initialization, fast)

    # loop until len(remain) == 0 or len(lVOP) >= max_vop_nb
    while len(remain) > 0 and len(lVOP) < max_vop_nb:
        # Compute sar criterion on the remaining Qs
        _rqstar_result = None
        if len(lVOP) > 0:
            if pool is None:
                _rqstar_result = [
                    _compute_rqstar(lQ[_remain],
                                    lVOP,
                                    real,
                                    max_spectral_norm,
                                    rqstar[1][_remain],
                                    rqstar_margin == 0)
                    for _remain in remain]
            else:
                # use parallel computing
                _rqstar_result = pool.starmap(
                    _compute_rqstar,
                    [(lQ[_remain],
                      lVOP,
                      real,
                      max_spectral_norm,
                      rqstar[1][_remain],
                      rqstar_margin == 0) for _remain in remain])
        # determine next vop candidate
        _max_spectral_norm = -np.Inf
        _max_rqstar = -np.Inf
        _vop_candidate = None
        for _i in range(len(remain)):
            _remain = remain[_i]
            if _rqstar_result is not None:
                # _max_rqstar
                rqstar[0][_remain] = _rqstar_result[_i][0]
                rqstar[1][_remain] = _rqstar_result[_i][1]
                # update VOP classif
                if rqstar[0][_remain] <= 1.0 + rqstar_margin:
                    classif[_remain] = 0
                # compute _max_rqstar and update _vop_candidate if condition
                # below is fulfilled
                if rqstar[0][_remain] > _max_rqstar:
                    _max_rqstar = rqstar[0][_remain]
                    if classif[_remain] < 0 and rqstar_margin > 0:
                        _vop_candidate = _remain
            # _max_spectral_norm
            if classif[_remain] < 0 and (
                    spectral_norm[_remain] > _max_spectral_norm):
                _max_spectral_norm = spectral_norm[_remain]
                _vop_candidate = _remain
        if _vop_candidate is None:
            print("No VOP to add anymore and _max_rqstar = %f" % _max_rqstar)
            break
        print(
            "%d VOPs / %d Q-matrices and %d to test still, "
            "_max_rqstar = %f, _max_spectr_norm = %f, new VOP candidate = %d"
            % (len(lVOP), N, len(remain), _max_rqstar, _max_spectral_norm,
               _vop_candidate))
        # increase VOP list and update classif
        lVOP.append(lQ[_vop_candidate] + compression_margin)
        classif[_vop_candidate] = 1
        # update 'remain' list
        remain = np.where(classif < 0)[0]
    if rqstar_margin > 0:
        for i in range(len(lVOP)):
            lVOP[i] *= _max_rqstar
    return lVOP, classif, rqstar


def _hermitian_to_realsym(Q):
    if isinstance(Q, list):
        return [_hermitian_to_realsym(_Q) for _Q in Q]
    else:
        return np.concatenate(
            (
                np.concatenate((np.real(Q), -np.imag(Q)), axis=-1),
                np.concatenate((np.imag(Q), np.real(Q)), axis=-1),
            ),
            axis=-2,
        )


def _compute_rqstar(Q, lVOP, real, scaling, initialization, fast):
    return compute_rqstar(
        Q,
        lVOP,
        real=real,
        scaling=scaling,
        initialization=initialization,
        fast=fast)


def compute_rqstar(Q,
                   lVOP,
                   real=None,
                   scaling=1.0,
                   initialization=None,
                   fast=False,
                   ):

    # Test if Q <= Qvop
    # input :
    #    Q10g : NcxNc Hermitian matric
    #    Qvop : list of Nvop NcxNc Hermitian matrices
    # Output :
    #    R       N real>0      max ((V^H Q V)/(V^H Qvop V))
    #    V       NcxN complex    argmax ((V^H Q V)/(V^H Qvop V))
    Nvop = len(lVOP)
    assert Q.ndim == 2, "Q is expected to be a single SAR matrix"
    Nc = Q.shape[0]
    for Qvop in lVOP:
        assert (
            Qvop.ndim == 2 and Qvop.shape[0] == Nc and Qvop.shape[1] == Nc
        ), "Bad list of VOP (shape)"
    # real symmetric ?
    if real is None:
        real = _is_real(Q)
    # Convert to real symetric matrices ?
    if real is False:
        lVOP = _hermitian_to_realsym(lVOP)
        Q = _hermitian_to_realsym(Q)
        Nc = 2 * Nc
        if initialization is not None:
            initialization = np.concatenate(
                (initialization.real, initialization.imag), axis=0)
    # Initialization vector
    if initialization is None:
        initialization = np.zeros((Nc,))
    if np.all(initialization == 0):
        initialization[0] = 1
        _, eigenvector = np.linalg.eigh(Q / scaling)
        initialization = eigenvector[:, Nc - 1]
        C = max(
            [
                _sar(Qvop, initialization, real=True, scaling=scaling)
                for Qvop in lVOP
            ]
        )
        initialization *= 0.5 / (np.sqrt(max(0, C)) + 1e-9)
    nlconstraints = scipy.optimize.NonlinearConstraint(
        lambda V: _constrfun_scipy(lVOP, V, scaling),
        np.zeros((Nvop)) - np.inf,
        np.zeros((Nvop)),
        jac=lambda V: _constrfun_jac_scipy(lVOP, V, scaling),
    )
    callback = None
    if fast:
        callback = lambda V: _stop_criterion_scipy(Q, lVOP, V, scaling)
    res = scipy.optimize.minimize(
        lambda V: _objfun_scipy(Q, V, scaling),
        initialization,
        args=(),
        method="SLSQP",
        jac=True,
        constraints=nlconstraints,
        callback=callback,
        tol=1e-8,
        options=dict(maxiter=1000, disp=False),
    )
    rqstar_maximizer = res.x
# Store solution of the problem  max ((V^H Q V)/(V^H Qvop V))
# rqstar_value = sar(Q, rqstar_maximizer, real=True, scaling=scaling) / max(
#    [sar(Qvop, rqstar_maximizer, real=True, scaling=scaling) for Qvop in lVOP]
# )
    rqstar_value = _sar(Q, rqstar_maximizer, real=True, scaling=scaling) / max(
        [_sar(Qvop, rqstar_maximizer, real=True, scaling=scaling)
         for Qvop in lVOP])
    if real is False:
        rqstar_maximizer = (
            rqstar_maximizer[0:Nc // 2] + 1j * rqstar_maximizer[Nc // 2:Nc]
        )
    return (rqstar_value, rqstar_maximizer)


def _sar(Qvop, V, real=None, scaling=1.):
    if real is None and (
        np.iscomplexobj(Qvop) or np.iscomplexobj(V)
    ):
        return V.T.conj() @ Qvop @ V / scaling
    else:
        return V.T @ Qvop @ V / scaling


# Objective function
def _objfun_nlopt(Q, V, G):
    P = -V.transpose() @ Q
    C = P @ V
    G[:] = 2 * P
    return C


def _objfun_scipy(Q, V, scaling):
    G = -V.transpose() @ Q / scaling
    C = G @ V
    G *= 2
    return (C, G)


# Constraints


def _constrfun(lQ, C, V, G):
    Vt = V.transpose()
    for i, Q in zip(range(0, len(lQ)), lQ):
        P = Vt @ lQ[i]
        C[i] = P @ V - 1
        G[i, :] = 2 * P


def _constrfun_scipy(lQ, V, scaling):
    Vt = V.transpose()
    C = np.zeros((len(lQ)))
    for i, Q in zip(range(0, len(lQ)), lQ):
        C[i] = Vt @ lQ[i] @ V / scaling - 1
    return C


def _constrfun_jac_scipy(lQ, V, scaling):
    Vt = V.transpose()
    G = np.zeros((len(lQ), V.size))
    for i, Q in zip(range(0, len(lQ)), lQ):
        G[i, :] = 2 * Vt @ lQ[i] / scaling
    return G


def _stop_criterion_scipy(Q, lVOP, V, scaling):
    fval = _objfun_scipy(Q, V, scaling)[0]
    constr = _constrfun_scipy(lVOP, V, scaling)
    constr = np.max(constr)
    return fval < -1 and constr <= 0


def _is_real(Q):
    if isinstance(Q, list):
        if any([_is_real(_Q) for _Q in Q]):
            assert all(
                [_is_real(_Q) for _Q in Q]
            ), "expecting either all real matrices or all complex"
            return True
        else:
            return False
    else:
        assert isinstance(Q, np.ndarray), "expecting nd-array"
        return not np.iscomplexobj(Q)


def _is_real_symmetric(Q, rtol=1e-5, atol=1e-8):
    if isinstance(Q, list):
        if any([_is_real_symmetric(_Q, rtol, atol) for _Q in Q]):
            assert all(
                [_is_real_symmetric(_Q, rtol, atol) for _Q in Q]
            ), "expecting either all real symmetric matrices or all Hermitian"
            return True
        else:
            return False
    else:
        assert isinstance(Q, np.ndarray), "expecting nd-array"
        Qt = Q.swapaxes(-1, -2)
        if np.iscomplexobj(Q):
            assert np.allclose(Q, Qt.conj(), rtol=rtol, atol=atol), (
                "Found non-Hermitian matrices")
            return False
        else:
            assert np.allclose(Q, Qt.conj(), rtol=rtol, atol=atol), (
                "Found non-symmetric matrices")
            return True
