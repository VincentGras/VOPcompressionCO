load testQ.mat Q10g

S = spectralNorm(Q10g);
Smax = max(S)

Qmargin = 0.05 * Smax * eye(8);

%% CO
    
[c, Qvop1] =  computeVOP_CO(Q10g, [], [], Qmargin, [], []);
check = rQstar(Q10g, Qvop1);
max(check)

%% iCO
niter = 10;
r = 0.8; 
[c, Qvop2, Qmargins] =  computeVOPi_CO(Q10g, Qmargin / (r^(niter-1)), r, niter);
Qvop2 = Qvop2{end};
check = rQstar(Q10g, Qvop2);
max(check)





