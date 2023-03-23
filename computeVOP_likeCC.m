function [c, Qvop] = computeVOP_likeCC(metric, Q, Qmargin, c, QvopBase)

% Compute a set of VOPs using the CO method
% Input:
%   criterion   function handle (ex: @(R,J) R(J))
%   Q       NcxNcxN complex Hermitian>0 ( or [] )
%   Qmargin NcxNcx1 complex Hermitian>0 (or positive scalar) : the
%              compression margin
%   c       1xN integer ( or [] )  initial classification
% Output :
%   c       1xN integer : the classification of each Q matrix, c(i) = 2 if VOP, 1 otherwise
%   Qvop    Updated set of VOPs, if Qvop is provided and Q <= Qvop (combined domination), then Qvop left unchanged
%           Note that any added VOP will be of the form Q(:,:,i) + Qmargin
%           for some i in {1,...,N}

% nb of channels
Nc = size(Q, 1);

% number od SAR matrices
N = size(Q, 3);

% External set of VOPs

if (nargin < 5)
    QvopBase = zeros(Nc,Nc,0);
end

% sort matrices

S = spectralNorm(Q);
[S, p] = sort(S, 'descend');
Q = Q(:,:,p);

% classif

[nyc, vop, dominated] = computeVOP_classif_code();

if (nargin < 4 || isempty(c))
    % initialize c as nyc (not yet classified)
    c = zeros(1,N);
    c(:) = nyc;
end

assert(ndims(c) <= 2 && numel(c) == N, 'Bad input c');

% VOP compression margin 

if (nargin < 3 || isempty(Qmargin))
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
    R(J) = metric(Q(:,:,J), cat(3, Qvop, QvopBase));
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

