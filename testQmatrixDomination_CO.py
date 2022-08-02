# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:14:52 2022

@author: VG244206
"""
import numpy as np
import nlopt

def testQmatrixDomination_CO(Q, lVOP):

    # Test if Q <= Qvop (combined domination)
    # input :
    #    Q10g : NcxNc Hermitian matric 
    #    Qvop : list of Nvop NcxNc Hermitian matrices
    # Output :
    #    R       N real>0      max ((V^H Q V)/(V^H Qvop V))
    #    V       NcxN complex    argmax ((V^H Q V)/(V^H Qvop V))

    
    Nvop = len(lVOP);            
    Nc = Q.shape[0];
    
    for Qvop in lVOP:
        assert Qvop.ndim == 2 and Qvop.shape[0]==Nc and Qvop.shape[1]==Nc, 'Bad list of VOP (shape)'
    
    
    # Convert to real symetric matrices
    CHtoRS = lambda Q: np.vstack((np.hstack((np.real(Q), -np.imag(Q))), np.hstack((np.imag(Q), np.real(Q)))));
    
    QvopR = [ CHtoRS(Qvop) for Qvop in lVOP];
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
    G[:] = 2 * P  #.flatten();
    
    return C;


# Constraints

def constrfun(lQ,C,V,G) :
    
    Vt = V.transpose();
    
    for i,Q in zip(range(0,len(lQ)),lQ) :
        
        P = Vt @ lQ[i];
        C[i] = P @ V - 1;
        G[i,:] = 2*P;

