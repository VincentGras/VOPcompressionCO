function [c, Qvop, Qmargin_i, ct] = computeVOPi(metric, Q, QmarginInitial, niter, r, max_vop_nb, QvopBase)

% VOP computation: iterative CO approach


Nc = size(Q, 1);
S = spectralNorm(Q);

% External set of VOPs

if (nargin < 7)
    QvopBase = zeros(Nc,Nc,0);
end

% max nb of VOPs

if (nargin < 6)
    max_vop_nb = inf;
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
    % set default initial margin to 0.5 max_Q(||Q||_2)
    QmarginInitial = max(S) * 0.5 * eye(Nc);
end

if (isscalar(QmarginInitial))
    QmarginInitial = QmarginInitial * eye(Nc);
end

% initialize VOP compression margin
Qmargin_i = QmarginInitial;

% initialize VOP classif
[nyc, vop, dominated] = computeVOP_classif_code();
c_i = zeros(1, size(Q, 3)) + nyc;

% prepare storage for VOP classif and VOP list at each iteration
c = zeros(niter, size(Q, 3));
Qvop = cell(niter, 1);
Qmargin = cell(niter, 1);

% ct

ct = zeros(niter, 1);

tic;

% iterate
for i = 1:niter
    % update classification
    Qmargin{i} = Qmargin;
    [c_i, Qvop{i}] = computeVOP(metric, Q, Qmargin_i, c_i, QvopBase);
    ct(i) = toc;
    c(i,:) = c_i;
    fprintf('Nb of VOP after iteration %d/%d = %d; (margin = %.1e)\n', i, niter, nnz(c_i), norm(Qmargin_i));
    
    if (nnz(c_i==vop) > max_vop_nb)
        fprintf('Stopping at iteration #%d as te nb of vops exceeds %d !', i, max_vop_nb);
        break;
    end
    
    % Decrease VOP compression margin
    Qmargin_i = Qmargin_i * r;
    
    % Update c
    c_i(c_i==dominated)=nyc;
end
