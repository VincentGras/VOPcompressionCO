function [c, Qvop] = computeVOP_SCO(Q, varargin)

% Compute a set of VOPs using the CO method and the criterion @(R,S) R
% See computeVOP_General_CO for details


[c, Qvop] = computeVOP_General_CO(@testQmatrixDomination_SCO, Q, spectralNorm(Q), varargin{:});
