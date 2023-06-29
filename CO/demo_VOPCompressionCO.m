load testQ.mat Q10g

S = spectralNorm(Q10g);
Smax = max(S)

Qmargin = 0.05 * Smax * eye(8);

%% CO
    
[c1, Qvop1] =  computeVOP_CO(Q10g, [], [], Qmargin, [], []);

% Equivalently we coult do (VOP not sorted according to spectral norm though):
% Qvop1 = Q10g(:, :, c1 > 0) + Qmargin;

%% iCO
niter = 10;
r = 0.8; 
[c2, Qvop2, Qmargins] =  computeVOPi_CO(Q10g, Qmargin / (r^(niter-1)), r, niter);
% here c2 is of size Niter x Nmatrices and Qvop2 is a cell array of size Niter !
% Fetch the result at last iteration :
Qvop2 = Qvop2{end};

% Equivalently we coult do (VOP not sorted according to spectral norm though):
% Qvop2 = Q10g(:, :, c2(end, :) > 0) + Qmargin(:, :, end);

%% Checks

check1 = rQstar(Q10g, Qvop1);
fprintf('N1* = %d, max(rQstar(Q10g, Qvop1)) = %f\n', sum(c1), max(check1));
check2 = rQstar(Q10g, Qvop2);
fprintf('N2* = %d, max(rQstar(Q10g, Qvop2)) = %f\n', sum(c2(end, :)), max(check2));
