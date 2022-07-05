function [Qvop, c] = computeVOP_CO(Q, varargin)

% Compute a set of VOPs using the CO method and the criterion @(R,S) R
% See computeVOP_General_CO for details


[Qvop, c] = computeVOP_General_CO(@(R,S) R, Q, spectralNorm(Q), varargin{:});
