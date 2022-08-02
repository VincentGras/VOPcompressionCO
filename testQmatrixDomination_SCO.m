function [R,V] = testQmatrixDomination_SCO(Q, Qvop)

% Test whether Q <= Qvop (combined domination)
% Input : 
%   Q       NcxNcxN complex Hermitian>0
%   Qvop    NcxNcxNvop complex Hermitian>0
% Output :
%   R       1xN real>0      max ((V^H Q V)/(V^H Qvop V))
%   V       NcxN complex    argmax ((V^H Q V)/(V^H Qvop V))
% Requirement :
%   MATLAB R2020b 
%   Optimization toolbox

maxiter = 1000;

N = size(Q, 3);
Nc = size(Qvop, 1);

[R,VCO] = testQmatrixDomination_CO(Q, Qvop);
VRCO = cat(1, real(VCO), imag(VCO));
V = zeros(Nc, Nc, N);
V(:,1,:) = reshape(VCO, Nc,1,N);

if (size(Qvop, 3) <= 1)
    return;
end

CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
%RStoCH = @(Q, Nc) Q(1:Nc,1:Nc) + 1i*Q(1+Nc:end,1:Nc);
QvopR = CHtoRS(Qvop);

[I,J] = indexing(2*Nc);
%[A,B] = ineqlinconstr(J) ;

opt = optimoptions(@fmincon, 'Algorithm', 'sqp', 'MaxIter', maxiter,  'Display', 'off', 'GradObj', 'on', 'GradConstr', 'on'); %,'OutputFcn',outfun);
if false
    opt = optimoptions(@fmincon, 'Algorithm', 'sqp', 'MaxIter', maxiter,  'Display', 'iter', 'GradObj', 'on', 'GradConstr', 'on', ...
        'CheckGradients',false);
end

parfor i = 1:N
    
    QtR = CHtoRS(Q(:,:,i));
    QtR = 0.5*(QtR+QtR');
        
%     [eigV, eigD] = eig(QtR, 'vector');
%     eigD = real(eigD);
%     [~, k] = max(eigD);
%     V0R = real(eigV(:,k)); % make sure real
%     %V0R(1) = 1; V0R(2:end)=0;
%     V0R = 0.5*V0R / (sqrt(max(evalSARReal(QvopR, V0R)))+eps);
    
%    U0R = zeros(size(I));
        %   U0R(1:2*Nc) = VRCO(:,i)/(sqrt(max(evalSARReal(QvopR, VRCO(:,i)))*0.95+eps));

    V0R = zeros(2*Nc);
    V0R(:,1) = VRCO(:,i);
    V0R = V0R + (eps + 1e-3 * norm(V0R(:,1))) + eye(2*Nc);
    V0R = V0R / sqrt(max(evalSARReal(QvopR, V0R))) ;
    U0R = V0R(I);
 %   max(constr(QvopR, U0R, I, J))
 %   max(A*U0R-B)
    
    [UR, ~, ~, output] = fmincon(@(UR) objfun(QtR,UR, I, J), U0R, [],[], [], [], [], [], @(UR) constr(QvopR, UR, I, J), opt);
    
    VR = zeros(2*Nc, 2*Nc);
    VR(I) = UR;
    
    oldR = R(i);
    newR = evalSARReal(QtR,VR)./max(evalSARReal(QvopR,VR));
    
    if (newR >= R(i))
        R(i) = newR;
        V(:,:,i) = VR(1:Nc,1:Nc) + 1i *  VR(1+Nc:end,1:Nc);
    end
    
    if (mod(i,100)== 0)
        fprintf('parfor loop : %d / %d; R/R-oldR= %f/%f; iterations = %d;\n', i, N, newR, newR-oldR,  output.iterations);
    end
    
    
end

function [I,J] = indexing(n)

    nn = n*n;
    m = (n+1)*n/2;
    I = zeros(m,1);
    J = zeros(n,n);
    L = reshape(1:nn,n,n);
    
    k=1;
    for j = 1:n
        for i = j:n
            I(k) = L(i,j);
            J(i,j) = k;
            k=k+1;
        end
    end



function [C, G] = objfun(Q, U, I, J)

n = size(J, 1);
nn=n*n;
V = zeros(n);
V(I) = U;
W = V*V';
V = sparse(V);

C = -trace(W * Q);
%G_ = zeros(size(U));
G = zeros(size(U));

% deriver V par rapport aux variables

for k = 1:numel(U)

    i = mod(I(k)-1, n)+1;
    j = 1+(I(k)-i)/n;
    G(k) = -2 * Q(i,:)*V(:,j);

%     nablaX1 = zeros(n);
%     nablaX1(i,:) = V(:,j)';
%     nablaX2 = zeros(n);
%     nablaX2(:,i) = V(:,j)';
    %nablaX_ = nablaX1+nablaX2;
    
        % nablaX = Colj' en i + Colj en i

    %for p = 1:n
        %G(k) = G(k) - nablaX1(i,p) * Q(p,i) - nablaX2(p,i) * Q(i, p);
        %G(k) = G(k) - V(p,j) * Q(p,i) - V(p,j) * Q(i, p);
        %G_(k) = G_(k) - 2 * V(p,j) * Q(p,i);
    %end

    continue;
    %[i,j] = ind2sub([n,n], I(k));
    nablaV = sparse(i,j,1,n,n);
    nablaX = nablaV*V'+V*nablaV';
    [ii,jj,nablaXval]=find(nablaX);
    
    
    
    for p = 1:numel(ii)
        G(k) = G(k) - nablaXval(p) * Q(jj(p), ii(p));
    end
    
end


% function [A,B] = ineqlinconstr(J) 
% 
% n = size(J, 1);
% A = zeros(n, n*(n+1)/2);
% B = zeros(n,1);
% for i = 1:n
%     A(i,J(i,i))=-1;
% end


function [C,Ceq, G, Geq] = constr(Q,U,I, J)

n = size(J, 1);
nn = n*n;

Ceq = [];
Geq = [];

V = zeros(n);
V(I) = U;
W = V*V';
C = zeros(size(Q,3),1);

for ic = 1:numel(C)
    C(ic) = trace(W*Q(:,:,ic))-1;
end

G = zeros(numel(U),size(Q,3));


for k = 1:numel(U)

    i = mod(I(k)-1, n)+1;
    j = 1+(I(k)-i)/n;
    
    for l = 1:size(Q,3)
        G(k,l) = 2 * Q(i,:,l)*V(:,j);
    end
% 
%     
%     [i,j] = ind2sub([n,n], I(k));
%     nablaV = sparse(i,j,1,n,n);
%     nablaX = nablaV*V'+V*nablaV';
%     [i,j,nablaXval]=find(nablaX);
%     
%     for p = 1:numel(i)
%         G(k,:) = G(k,:) + nablaXval(p) * reshape(Q(j(p), i(p), :), 1, size(Q,3));
%     end
    
end



function SED = evalSARReal(Q, V)

SED = zeros(size(Q,3),1);
W = V*V';
for i = 1:numel(SED)
    SED(i) = max(trace(W*Q(:,:,i)),0); 
end




