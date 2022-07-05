function [Qvop,c] = computeVOP_iCO(criterion, Q, Qmargin, niter, r)

% VOP computation: iterative CO approach
% Input :
%   criterion       function handle @(R,S) R or @(R,S) S
%   Q               NcxNcxN complex Hermitian>0
%   margin          scalar > 0 : desired margin (default = 0.01 max_Q(||Q||_2))
%   niter           integer    : number of iterations (default = 4)
%   r               scalar (0<r<1) : margin update factor (default = 0.6)
% Output :
% same as computeVOP_CO

Nc = size(Q, 1);
S = spectralNorm(Q);

if (nargin < 5)
    % default r
    r = 0.8;
end

assert(r>0 && r < 1, 'bad input parameter r, must be with the range ]0,1[');

if (nargin < 4)
    niter = 4;
end

if (nargin < 3)
    % set default margin to 0.01 max_Q(||Q||_2)
    Qmargin = max(S) * 0.05 * eye(Nc);
end

if (isscalar(Qmargin))
    Qmargin = Qmargin * eye(Nc);
end

% intialize Qvop to []
Qvop = [];

% initialize VOP compression margin
Qm = Qmargin / (r^(niter-1));
fprintf('Initial margin = %.1e\n', norm(Qm));

% initialize c
[nyc, ~, dominated] = computeVOP_classif_code();
c = zeros(1, size(Q, 3));
c(:) = nyc;


% iterate
for i = 1:niter
    % update classification
    [Qvop, c] = computeVOP_General_CO(criterion, Q, S, Qm, c);
    fprintf('Nb of VOP after iteration %d/%d = %d; (margin = %.1e)\n', i, niter, nnz(c), norm(Qm));
    % Decrease VOP compression margin
    Qm = Qm * r;
    % Update c
    c(c==dominated)=nyc;
end


function [SARwc] = maxSpectralNorm(Q)

N = size(Q, 3); 

SARwc = zeros(N,1);

parfor i = 1:N
    SARwc(i) = norm(Q(:,:,i), 2);
end

SARwc = max(SARwc);
