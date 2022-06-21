function [Qvop,c] = computeVOP_iCO(Q, Qmargin, niter, r)

% VOP computation: iterative CO approach
% Input :
%   Q               NcxNcxN complex Hermitian>0
%   margin          scalar > 0 : desired margin (default = 0.01 max_Q(||Q||_2))
%   niter           integer    : number of iterations (default = 4)
%   r               scalar (0<r<1) : margin update factor (default = 0.6)
% Output :
% same as computeVOP_CO

Nc = size(Q, 1);
wcSAR = maxSpectralNorm(Q);

if (nargin < 4)
    % default r
    r = 0.6;
end

assert(r>0 && r < 1, 'bad input parameter r, must be with the range ]0,1[');

if (nargin < 3)
    niter = 4;
end

if (nargin < 2)
    % set default margin to 0.01 max_Q(||Q||_2)
    Qmargin = wcSAR * 0.01 * eye(Nc);
end

if (isscalar(Qmargin))
    Qmargin = Qmargin * eye(Nc);
end

% intialize Qvop to []
Qvop = [];

% initialize VOP compression margin
Qm = Qmargin / (r^(niter-1));

% initialize c
c = zeros(1, size(Q, 3));

% iterate
for i = 1:niter
    % update classification
    [Qvop, c] = computeVOP_CO(Q, Qm, Q(:, :, c > 0)+Qm);
    fprintf('Nb of VOP after iteration %d/%d = %d; (margin = %.1e)\n', i, niter, nnz(c), norm(Qm));
    % Decrease VOP compression margin
    Qm = Qm * r;
end


function [SARwc] = maxSpectralNorm(Q)

N = size(Q, 3); 

SARwc = zeros(N,1);

parfor i = 1:N
    SARwc(i) = norm(Q(:,:,i), 2);
end

SARwc = max(SARwc);
