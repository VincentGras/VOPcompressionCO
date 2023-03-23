function [c, Qvop, Qmargin_i, ct] = computeVOPi_CO(Q, QmarginInitial, niter, r, max_vop_nb, QvopBase)

% VOP computation: iterative CO approach

Nc = size(Q, 1);

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

% nb of iterations

if (nargin < 4)
    niter = 4;
end

% initial Qmargin

if (nargin < 3)
    QmarginInitial = [];
end


if (~isreal(Q))
    assert(isreal(QmarginInitial), 'Q is real but not QmarginInitial');
    assert(isreal(QvopBase), 'Q is real but not QvopBase');
    % transform to symmetric Hermitian and recursive call
    CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
    Q = CHtoRS(Q);
    if ~isempty(QmarginInitial) && ~isscalar(QmarginInitial)
        QmarginInitial = CHtoRS(QmarginInitial);
    end
    if ~isempty(QvopBase)
        QvopBase = CHtoRS(QvopBase);
    end
    [c, Qvop, Qmargin_i, ct] = computeVOPi_CO(Q, QmarginInitial, niter, r, max_vop_nb, QvopBase);
    for i = 1:numel(Qvop)
        Qvop{i} = Qvop{i}(1:Nc, 1:Nc, :) + 1i * Qvop{i}(Nc+1:2*Nc, 1:Nc, :);
    end
    return;
end

[S, p] = sort(spectralNorm(Q), 'descend');
Q = Q(:,:,p);

if (isempty(QmarginInitial))
    % set default initial margin to 0.5 max_Q(||Q||_2)
    QmarginInitial = S(1) * 0.5 * eye(Nc);
end

if isscalar(QmarginInitial)
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

X=[]; 

% iterate
for i = 1:niter
    % update classification
    Qmargin{i} = Qmargin;
    [c_i, Qvop{i}, X] = computeVOP_CO(Q, S, X, Qmargin_i, c_i, QvopBase);
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

% inverse permutation
c_old = c;
invp = 1:numel(p); 
invp(p)=invp;
c = c(:, invp);
figure; plot(c)

