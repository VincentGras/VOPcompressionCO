function [Qvop,c] = computeVOP_iCO(criterion, Q, varargin{:})

% VOP computation: iterative CO approach with criterion = R
% See also  computeVOP_General_iCO

[Qvop,c] = computeVOP_General_iCO(@(R,J) R(J), Q, varargin);