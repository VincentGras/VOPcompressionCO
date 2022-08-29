function [M, x] = computeMarginSAR(Q, Qvop)


[M, X] = computeCriterionSAR(Q, Qvop);
M = 1./M - 1;
