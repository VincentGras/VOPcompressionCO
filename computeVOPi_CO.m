function [c, Qvop, Qmargin] = computeVOPi_CO(Q, QmarginInitial, r, max_iter, max_vop, QvopInitial)

% VOP computation: iterative CO approach
% Q :               Nc x Nc x N complex or real
% QmarginInitial :  Nc x Nc or scalar
%                   initial compression margin
% r              :  iteration factor (0 < r < 1)
% max_iter       : max nb of iteration
% max_vop        : max vop number (default Inf)
% QvopInitial    :  Initial VOP set 
%                   default is []

Nc = size(Q, 1);

% External set of VOPs

if (nargin < 6)
    QvopInitial = zeros(Nc,Nc,0);
end

if (nargin < 5)
    max_vop = inf;
end

assert(r>0 && r < 1, 'bad input parameter r, must be with the range ]0,1[');

if (~isreal(Q))
    % Convert complex Hermitian matrices into symmetric real matrices
    % transform to symmetric Hermitian and recursive call
    CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
    Q = CHtoRS(Q);
    if ~isempty(QmarginInitial) && ~isscalar(QmarginInitial)
        QmarginInitial = CHtoRS(QmarginInitial);
    end
    if ~isempty(QvopInitial)
        QvopInitial = CHtoRS(QvopInitial);
    end
    [c, Qvop, Qmargin] = computeVOPi_CO(Q, QmarginInitial, r, max_iter, max_vop, QvopInitial);
    Qmargin = Qmargin(1:Nc, 1:Nc, :) + 1i * Qmargin(Nc+1:2*Nc, 1:Nc, :);
    Qvop = cellfun(@(Q) Q(1:Nc, 1:Nc, :) + 1i * Q(Nc+1:2*Nc, 1:Nc, :), Qvop, 'UniformOutput', false);
    return;
end

% check that all matrices are real
assert(isreal(QmarginInitial), 'Q is real but not QmarginInitial');
assert(isreal(QvopInitial), 'Q is real but not QvopBase');

% compute spectral norm across all Qs and sort
% keep permutation p for later
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
c = zeros(max_iter, size(Q, 3));
Qvop = cell(max_iter, 1);
Qmargin = zeros(Nc, Nc, max_iter);

% RF shim maximizing xHQx subbject to VOP constraints
X=[]; 

% computation time
%ct = zeros(max_iter, 1);
%tic;

% iterate
for iter = 1:max_iter
    % update classification
    Qmargin(:, :, iter) = Qmargin_i;
    [c_i, Qvop{iter}, X] = computeVOP_CO(Q, S, X, Qmargin_i, c_i, QvopInitial);
    %ct(iter) = toc;
    c(iter,:) = c_i;
    fprintf('Nb of VOP after iteration %d/%d = %d; (margin = %.1e)\n', iter, max_iter, nnz(c_i), norm(Qmargin_i));
    
    if (nnz(c_i==vop) > max_vop)
        fprintf('Stopping at iteration #%d as te nb of vops exceeds %d !', iter, max_vop);
        break;
    end
    
    % Decrease VOP compression margin
    Qmargin_i = Qmargin_i * r;
    
    % Update c
    c_i(c_i==dominated)=nyc;
end

if (iter < max_iter)
    c = c(1:iter, :);
    Qmargin = Qmargin(:, :, 1:iter);
    Qvop = Qvop(1:iter);
    %ct = ct(1:niter);
end

% inverse permutation
invp = 1:numel(p); 
invp(p)=invp;
c = c(:, invp);

