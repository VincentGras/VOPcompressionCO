function [c, Qvop] = computeVOP_iSCO(Q, varargin)

% VOP computation: iterative CO approach with criterion = R
% See also  computeVOP_General_iCO

[c, Qvop] = computeVOP_General_iCO(@testQmatrixDomination_SCO, Q, varargin{:});