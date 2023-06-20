function [c, Qvop, X] = computeVOP_CO(Q, S, X, Qmargin, c, QvopInitial)

% Compute a set of VOPs using the CO method
% Input:
%   criterion   function handle (ex: @(R,J) R(J))
%   Q       NcxNcxN complex Hermitian>0 ( or [] )
%   S       Spectral norm of the Q matrices
%           If [], S is computed internally
%           !! If not [], it is assumed that S is sorted (and Q accordingly)
%   X       Nc x N array
%           First guess for the vector X maximazing the domination
%           criterion
%   Qmargin NcxNcx1 complex Hermitian>0 (or positive scalar) : the
%              compression margin
%   c       1xN integer ( or [] )  initial classification
% Output :
%   c       1xN integer : the classification of each Q matrix, c(i) = 2 if VOP, 1 otherwise
%   Qvop    Updated set of VOPs, if Qvop is provided and Q <= Qvop (combined domination), then Qvop left unchanged
%           Note that any added VOP will be of the form Q(:,:,i) + Qmargin
%           for some i in {1,...,N}
%   X       argmax of the domination criterion
%
%     CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
%     [CSAR, xr] = RQstar(CHtoRS(Q), [real(x); imag(x)], CHtoRS(Qvop), optimalityTolerance);
%     x = xr(1:Nc,:) + 1i * xr(Nc+1:2*Nc, :);
%
% else

% nb of channels
Nc = size(Q, 1);

if (~isreal(Q))
    assert(isreal(Qmargin), 'Q is real but not QmarginInitial');
    assert(isreal(QvopInitial), 'Q is real but not QvopBase');
    % transform to symmetric Hermitian and recursive call
    if (~isscalar(Qmargin))
        Qmargin = CHtoRS(Qmargin);
    end
    [c, Qvop, X] = computeVOP_CO(CHtoRS(Q), S, [real(X); imag(X)], Qmargin, c, QvopInitial);
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

if (isempty(QvopInitial))
    QvopInitial = zeros(Nc,Nc,0);
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

if (isempty(Qvop) && isempty(QvopInitial))
    % add first Q matrix (the one with the highest spectral norm)
    c(1) = vop;
    Qvop = Q(:,:,1) + Qmargin;
end

% loop

remain = nnz(c==nyc);

while (remain > 0)
    
    fprintf('%d VOPs / %d Q-matrices / %d to test still\n', nnz(c==vop), N, remain);
    
    J = find(c == nyc);
    
    % compute R on the nyc SAR matrices
	% Note the noramlization to Smax (i.e. S(1) as S is sorted) 
	% This is done to condition properly the numerical optimization
    [R(J), X(:, J)] = rQstar(Q(:,:,J)/S(1), cat(3, Qvop, QvopInitial)/S(1), true, X(:, J));
    
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

