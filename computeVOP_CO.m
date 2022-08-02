function [c, Qvop] = computeVOP_CO(Q, varargin)

% Compute a set of VOPs using the CO method and the criterion @(R,S) R
% See computeVOP_General_CO for details


[c, Qvop] = computeVOP_General_CO(@testQmatrixDomination_CO, Q, spectralNorm(Q), varargin{:});
