function [c, Qvop] = computeVOPi(metric, Q, Qmargin, niter, r, QvopExt)

% VOP computation: iterative CO approach


Nc = size(Q, 1);
S = spectralNorm(Q);

% External set of VOPs

if (nargin < 6)
    QvopExt = zeros(Nc,Nc,0);
end

% Qmargin update factor

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
    [c, Qvop] = computeVOP(metric, Q, S, Qm, c, QvopExt);
    fprintf('Nb of VOP after iteration %d/%d = %d; (margin = %.1e)\n', i, niter, nnz(c), norm(Qm));
    % Decrease VOP compression margin
    Qm = Qm * r;
    % Update c
    c(c==dominated)=nyc;
end
