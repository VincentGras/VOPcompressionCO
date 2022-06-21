# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:14:52 2022

@author: VG244206
"""
import numpy as np
import nlopt

def testQmatrixDomination_CO(Q, Qvop):

    # Test if Q <= Qvop (combined domination)
    # input :
    #    Q10g : NcxNxN ndarray   the 10g-SAR matrices
    #    Qvop : NcxNcxNvop VOPs 
    # Output :
    #    R       1xN real>0      max ((V^H Q V)/(V^H Qvop V))
    #    V       NcxN complex    argmax ((V^H Q V)/(V^H Qvop V))

    
    N = 1;
    if (Q.ndim>2) :
        N = Q.shape[2];
            
    Nc = Qvop.shape[0];
    
    V = np.zeros((Nc, N),dtype=complex);
    R = np.zeros(N);
    
    # complex 8*8 vers real 16x16 car nlopt ne travaille que sur des réels
    CHtoRS = lambda Q: np.vstack((np.hstack((np.real(Q), -np.imag(Q))), np.hstack((np.imag(Q), np.real(Q)))));
    
    Nvop = 1
    QvopR = CHtoRS(Qvop);
        
    if (QvopR.ndim>2):
        Nvop=QvopR.shape[2];
    
    if (Nvop == 1 and QvopR.ndim>2):
        QvopR = QvopR.reshape((2*Nc,2*Nc));    
    
    if Q.ndim<3:
        QR = CHtoRS(Q);
    
    for i in range(0,N) :
        
        if Q.ndim>2 :
            QR = CHtoRS(Q[:,:,i]);
        
        if Nvop == 1:
            
            # in this case do an eignevalue analysis of QR-QvopR
            # and look for the maximum eigenvalue
            # If that eigenvalue is negative then QR<=QvopR and R[i] will be < 1
             
            H = QR-QvopR;
            eigW,eigV = np.linalg.eigh(H);
            VR = eigV[:,2*Nc-1]; # vector corresponding to the max eigenvalue
            
        else:
            
            # In this case run a convex optimization
            
            V0R = np.zeros((2*Nc));
            V0R[0] = 1;
            C = np.max(evalSAR(QvopR, V0R));
            V0R *= (1.0 / np.sqrt(C)) ;
            
            # prepare nlopt     
            opt = nlopt.opt(nlopt.LD_MMA , 2*Nc)
            opt.set_min_objective(lambda x,grad: objfun(QR, x, grad))
            opt.add_inequality_mconstraint(lambda res,x,grad: VopConstr(QvopR, res, x, grad), 1e-8*np.ones(Nvop))
            opt.set_xtol_rel(1e-6)
            opt.set_maxeval(200)
            
            # find minimizer
            VR = opt.optimize(V0R)
        
        # Store solution of the problem  max ((V^H Q V)/(V^H Qvop V))
        R[i] = evalSAR(QR,VR) / evalSAR(QvopR,VR);
        V[:,i] = VR[0:Nc] + 1j * VR[Nc:2*Nc];

    return (R,V)

# Objective function

def objfun(Q, V, G) :
    P = -V.transpose() @ Q
    C = P @ V;
    G[:] = 2 * P  #.flatten();
    
    return C;


# Constraints

def VopConstr(Q,C,V,G) :
    
    Vt = V.transpose();
    
    if (Q.ndim<3):
        P = Vt @ Q;
        C[0] = P @ V - 1;
        G[:,0] = 2 * P;
    else :
        for i in range(0,Q.shape[2]) :
            
            P = Vt @ Q[:,:,i];
            C[i] = P @ V - 1;
            G[i,:] = 2*P;
    
# Evaluation of V^H Q V
 
def evalSAR(Q, V):
    
    Vt = V.transpose();
    
    if (Q.ndim<3):
        P = Vt @ Q;
        return P @ V;
        
    else :
        
        C = 0;
        
        for i in range(0,Q.shape[2]) :
            
            P = Vt @ Q[:,:,i];
            C = np.max((C, P @ V));
            
    return C