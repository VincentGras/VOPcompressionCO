function [c, Qvop, X] = computeVOP_CO(Q, S, X, Qmargin, c, QvopBase)

% Compute a set of VOPs using the CO method
% Input:
%   criterion   function handle (ex: @(R,J) R(J))
%   Q       NcxNcxN complex Hermitian>0 ( or [] )
%   S
%   X
%   Qmargin NcxNcx1 complex Hermitian>0 (or positive scalar) : the
%              compression margin
%   c       1xN integer ( or [] )  initial classification
% Output :
%   c       1xN integer : the classification of each Q matrix, c(i) = 2 if VOP, 1 otherwise
%   Qvop    Updated set of VOPs, if Qvop is provided and Q <= Qvop (combined domination), then Qvop left unchanged
%           Note that any added VOP will be of the form Q(:,:,i) + Qmargin
%           for some i in {1,...,N}

%
%     CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
%     [CSAR, xr] = RQstar(CHtoRS(Q), [real(x); imag(x)], CHtoRS(Qvop), optimalityTolerance);
%     x = xr(1:Nc,:) + 1i * xr(Nc+1:2*Nc, :);
%
% else

% nb of channels
Nc = size(Q, 1);

if (~isreal(Q))
    % transform to symmetric Hermitian and recursive call
    CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
    [c, Qvop, X] = computeVOP_CO(Q, S, [real(X); imag(X)], CHtoRS(Qmargin), c, QvopBase);
    X = X(1:Nc,:) + 1i * X(Nc+1:2*Nc, :);
    Qvop = Qvop(1:Nc, 1:Nc, :) + 1i * Qvop(Nc+1:2*Nc, 1:Nc, :);
    return;
end

% number od SAR matrices
N = size(Q, 3);

if (isempty(X))
    X = zeros(Nc, N);
end

% External set of VOPs

if (isempty(QvopBase))
    QvopBase = zeros(Nc,Nc,0);
end

% sort matrices ?

if (isempty(S))
    S = spectralNorm(Q);
    [S, p] = sort(S, 'descend');
    Q = Q(:,:,p);
end

% classif

[nyc, vop, dominated] = computeVOP_classif_code();

if (isempty(c))
    % initialize c as nyc (not yet classified)
    c = zeros(1,N);
    c(:) = nyc;
end

assert(ndims(c) <= 2 && numel(c) == N, 'Bad input c');

% VOP compression margin

if (isempty(Qmargin))
    %S = spectralNorm(Q);
    Qmargin = eye(Nc)*(0.05*max(S));
end

if (isscalar(Qmargin))
    Qmargin = Qmargin * eye(Nc);
end

% initialize R as inf
R = zeros(1,N) + inf;

% initialize Qvop

Qvop = Q(:,:,c==vop) + Qmargin;

if (isempty(Qvop) && isempty(QvopBase))
    Qvop = Qmargin;
end

% loop

remain = nnz(c==nyc);

while (remain > 0)
    
    fprintf('%d VOPs / %d Q-matrices / %d to test still\n', nnz(c==vop), N, remain);
    
    J = find(c == nyc);
    
    % compute R on the nyc SAR matrices
    [R(J), X(:, J)] = RQstar(Q(:,:,J), X(:, J), cat(3, Qvop, QvopBase));
    %R_ = testQmatrixDomination_CHO(Q(:,:,J), cat(3, Qvop, QvopExt));
    
    % mark as 'dominated' the matrices for which Rmax < 1
    c(J(R(J)<=1)) = dominated;
    J = find(c == nyc);
    
    % if  max(R)>1, mark as new vop the SAR matrix that realizes the max
    if (~isempty(J))
        %[~,k] = max(R(J));
        k = 1;
        c(J(k)) = vop;
        Qvop = Q(:,:,c==vop) + Qmargin;
        remain = numel(J) - 1;
    else
        remain = 0;
    end
    
end

fprintf('%d VOPs / %d Q-matrices \n', size(Qvop, 3), N);

function [CSAR, X] = RQstar(Q, X, Qvop)

Nc = size(Qvop, 1);
N = size(Q, 3);

assert(isreal(Q), 'Q must be real');
assert(isreal(Qvop), 'Qvop must be real');

% x = zeros(Nc, N);
CSAR = zeros(1, N);

parfor i = 1:N
    
    Qi = double(Q(:,:,i));
    Qi = 0.5*(Qi+Qi');
    
    X0 = X(:,i);
    
    if (all(X0 == 0))
        [eigV, eigD] = eig(Qi - Qvop(:, :, 1), 'vector');
        eigD = real(eigD);
        [~, k] = max(eigD);
        X0 = real(eigV(:,k));
    end
    
    X0 = sqrt( 1.0 / (SAR(Qvop, X0) + eps)) * X0;
    
    %         opt = optimset('Algorithm', 'interior-point', 'MaxIter', 1000, ...
    %             'Display', 'off', 'GradObj', 'on', 'GradConstr', 'on', ...
    %             'HessFcn', @(x,lambda) HessianFcn(Qi, Qvop, x, lambda), ...
    %             'OutputFcn', @stopCriterion);
    
    opt = optimset('Algorithm', 'sqp', 'MaxIter', 1000, ...
        'Display', 'off', 'GradObj', 'on', 'GradConstr', 'on', ...
        'OutputFcn', @stopCriterion);
    
    X(:,i) = fmincon(@(x) objfun(Qi,x), X0, [], [], [], [], [], [], @(x) constrfun(Qvop, x), opt);
    
    CSAR(i) = SAR(Qi,X(:,i))./SAR(Qvop,X(:,i));
    
    if (mod(i,1000)== 0)
        fprintf('parfor loop : %d / %d; %f\n', i, N, CSAR(i));
    end
    
    
end

function val = SAR(Q,V)
val = pagemtimes(V, 'transpose', pagemtimes(Q, V), 'none');
val = max(val, [], 'all');
val = max(val, 0);

function [C, G] = objfun(Q, V)
C = -V' * Q * V;
G = -2 * Q * V;

function [C,Ceq,G,Geq] = constrfun(Q,V)

Ceq = [];
C = reshape(pagemtimes(V, 'transpose', pagemtimes(Q, V), 'none'), size(Q, 3), 1) - 1;
Geq = [];
G = 2 * reshape(pagemtimes(Q, V), numel(V), size(Q, 3));

function H = HessianFcn(Q, Qvop, V, lambda)

H = -2 * Q;

for c = 1:size(Qvop, 3)
    H = H + 2 * lambda.ineqnonlin(c) * Qvop(:,:,c);
end

function stop = stopCriterion(x,optimValues,state)

stop = optimValues.fval < -1;

