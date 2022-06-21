function [Qvop,c] = computeVOP_CO(Q, Qmargin, Qvop)

% Compute a set of VOPs using the CO method
% Input:
%   Q       NcxNcxN complex Hermitian>0
%   Qmargin NcxNcx1 complex Hermitian>0 (or positive scalar) : the
%              compression margin
%   Qvop    NcxNcxNvop complex Hermitian>0 or []   [Optional]
% Output :
%   Qvop    Updated set of VOPs, if Qvop is provided and Q <= Qvop (combined domination), then Qvop left unchanged
%           Note that any added VOP will be of the form Q(:,:,i) + Qmargin
%           for some i in {1,...,N}
%   c       1xN integer : the classification of each Q matrix, c(i) = 2 if VOP, 1 otherwise

% nb of channels
Nc = size(Q, 1);

% number od SAR matrices
N = size(Q, 3);

% VOP compression margin 
if (nargin < 2)
    Qmargin = eye(Nc)*(0.01*maxSpectralNorm(Q));
end

if (isscalar(Qmargin))
    Qmargin = Qmargin * eye(Nc);
end

% Initial set of VOP (default = none)
if (nargin < 3)
    Qvop = [];
end

% values for c
dominated = 0;
vop = 1;
nyc = 2;

% initialize c as nyc (not yet classified)
c = nyc+zeros(1,N);

% initialize R as inf
R = zeros(1,N) + inf;

remain = N;

while (remain > 0)
    
    fprintf('%d new VOPs / %d Q-matrices / %d to test still\n', nnz(c==vop), N, remain);
  
    J = find(c == nyc);
    
    % compute R on the nyc SAR matrices
    if (isempty(Qvop))
        R(J) = testQmatrixDomination_CO(Q(:,:,J), Qmargin);
    else
        R(J) = testQmatrixDomination_CO(Q(:,:,J), Qvop);
    end
    
    % mark as 'dominated' the matrices for which Rmax < 1
    c(J(R(J)<=1)) = dominated;
    J = find(c == nyc);
    
    % if  max(R)>1, mark as new vop the SAR matrix that realizes the max 
    if (~isempty(J))
        [~,k] = max(R(J));
        c(J(k)) = vop;
        if (isempty(Qvop))
            Qvop = Q(:,:,J(k)) + Qmargin;
        else
            Qvop = cat(3, Q(:,:,J(k)) + Qmargin, Qvop);
        end
        remain = numel(J) - 1;
    else
        remain = 0;
    end
    
end

fprintf('%d new VOPs / %d Q-matrices \n', nnz(c==vop), N);


function [SARwc] = maxSpectralNorm(Q)

N = size(Q, 3); 

SARwc = zeros(N,1);

parfor i = 1:N
    SARwc(i) = norm(Q(:,:,i), 2);
end

SARwc = max(SARwc);

