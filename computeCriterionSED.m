function [CSED, X] = computeCriterionSED(Q, Qvop, optimalityTolerance)

if (nargin < 3)
    optimalityTolerance = 1e-6;
end

Nc = size(Qvop, 1);
N = size(Q, 3);

if (~isreal(Q) || ~isreal(Qvop))
    
    CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
    [CSED, xr] = computeCriterionSED(CHtoRS(Q), CHtoRS(Qvop), optimalityTolerance);
    X = xr(1:Nc,1:Nc,:) + 1i * xr(1+Nc:2*Nc, 1:Nc, :);
    
else
    
    X = zeros(Nc, Nc, N);
    CSED = zeros(1, N);
    
    [I,J] = indexing(Nc);
    
    
   parfor i = 1:N
        
        Qi = double(Q(:,:,i));
        Qi = 0.5*(Qi+Qi');
        
 
        if (size(Qvop, 3)==1)
            
            [eigV, eigD] = eig(Qi-Qvop, 'vector');
            eigD = real(eigD);
            [~, k] = max(eigD);
            Xi = zeros(Nc);
            Xi(:,1) = real(eigV(:,k));
            
        else
            Xi = chol(Qi + (1e-3*norm(Qi))*eye(Nc), 'lower');
            Xi =  sqrt(0.5 / (SED(Qvop, Xi)+eps))* Xi;   
            opt = optimoptions(@fmincon, 'Algorithm', 'interior-point', ...
                'MaxIter', 1000,  'Display', 'off', ...
                'GradObj', 'on', 'GradConstr', 'on', ...
                'OptimalityTolerance', optimalityTolerance, ...
                'StepTolerance', 1e-10, ...
                'FunctionTolerance', 1e-6, ...
                'HessianFcn', @(x,lambda) HessianFcn(Qi, Qvop, x, I, J, lambda));

            [Xi(I), fval, exitflag, output]  = fmincon(@(x) objfun(Qi,x, I, J), Xi(I), [], [], [], [], [], [], @(x) constrfun(Qvop, x, I, J), opt);
            
        end
        
        X(:,:,i) = Xi;
        CSED(i) = SED(Qi,X(:,:,i))/SED(Qvop,X(:,:,i));
        
        if (mod(i,1000) == 1)
            fprintf('parfor loop : %d / %d; %f\n', i, N, CSED(i));
        end
        
        
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
V = zeros(n);
V(I) = U;
W = V*V';
%V = sparse(V);

C = -trace(W * Q);
G = zeros(size(U));

% deriver V par rapport aux variables

% for k = 1:numel(U)
% 
%     i = mod(I(k)-1, n)+1;
%     j = 1+(I(k)-i)/n;
%     G(k) = -2 * Q(i,:)*V(:,j);
%     
% end

A = cumsum([0, n:-1:1]);

for j = 1:n
  
  G(A(j)+1:A(j+1)) = -2 * Q(j:n, :) * V(:, j);
    
end


function [C,Ceq, G, Geq] = constrfun(Q,U,I, J)

n = size(J, 1);

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

% li = mod(I-1, n)+1;
% co = 1+(I-li)/n;
% 
% for m = 1:size(Q,3)
%     for k = 1:numel(U)
%         G(k,m) = 2 *  Q(li(k),:,m)*V(:,co(k));
%     end
% end

A = cumsum([0, n:-1:1]);

% for m = 1:size(Q,3)
%     for j = 1:n
%         G(A(j)+1:A(j+1),m) = 2 * Q(j:n,:,m)*V(:,j);
%     end
% end

for j = 1:n
    G(A(j)+1:A(j+1),:) = 2 * pagemtimes(Q(j:n,:,:), V(:,j));
end


function H = HessianFcn(Q, Qvop, U, I, J, lambda)

p = numel(U);
n = size(J, 1);
H = zeros(p, p);

%Iinv = zeros(n*n, 1);
%Iinv(I) = 1:p;

A = cumsum([0, n:-1:1]);

% deriver V par rapport aux variables

lam = reshape(lambda.ineqnonlin, 1, 1, size(Qvop, 3));
Qvopxlam = sum(Qvop .* lam, 3);

for iH = 1:numel(U)

    iQ = mod(I(iH)-1, n) + 1;
    jQ = I(iH) - iQ;
    iA = 1 + jQ / n;
    jH = A(iA)+1:A(iA+1);
    h_ = -2 * Q(iQ, I(jH) - jQ) + 2 * Qvopxlam(iQ, I(jH) - jQ); 
    H(iH, jH) = h_;
end



function val = SED(Q, V)

val = 0;
W = V*V';
for i = 1:size(Q,3)
    val = max(trace(W*Q(:,:,i)),val); 
end

