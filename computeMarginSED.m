function [M, X] = computeMarginSED(Q, Qvop)


[M, X] = computeCriterionSED(Q, Qvop);
M = 1./M - 1;
