function [CSED, X] = computeCriterionSED(Q, Qvop)

Nc = size(Qvop, 1);
N = size(Q, 3);

if (~isreal(Q) || ~isreal(Qvop))
    
    CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
    [CSED, xr] = computeCriterionSED(CHtoRS(Q), CHtoRS(Qvop));
    X = xr(1:Nc,1:Nc,:) + 1i * xr(1:Nc, Nc+1:2*Nc, :);
    
else
    
    X = zeros(Nc, Nc, N);
    CSED = zeros(1, N);
    
    [I,J] = indexing(Nc);
    
    opt = optimset('Algorithm', 'sqp', 'MaxIter', 1000,  'Display', 'off', 'GradObj', 'on', 'GradConstr', 'on');
    
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
            Xi =  sqrt(0.5 / (SED(Qvop, x0)+eps))* Xi;
            Xi(I) = fmincon(@(x) objfun(Qi,x, I, J), Xi(I), [], [], [], [], [], [], @(x) constrfun(Qvop, x, I, J), opt);
            
        end
        
        X(:,:,i) = Xi;
        CSED(i) = SED(Qi,X(:,:,i))/SED(Qvop,X(:,:,i));
        
        if (mod(i,1000)== 0)
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
V = sparse(V);

C = -trace(W * Q);
G = zeros(size(U));

% deriver V par rapport aux variables

for k = 1:numel(U)

    i = mod(I(k)-1, n)+1;
    j = 1+(I(k)-i)/n;
    G(k) = -2 * Q(i,:)*V(:,j);
    
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

li = mod(I-1, n)+1;
co = 1+(I-li)/n;

for m = 1:size(Q,3)
    for k = 1:numel(U)
        G(k,m) = 2 *  Q(li(k),:,m)*V(:,co(k));
    end
end






function val = SED(Q, V)

val = 0;
W = V*V';
for i = 1:size(Q,3)
    val = max(trace(W*Q(:,:,i)),val); 
end


% function [C, G] = objfunMarginSAR(Q, V)
% C = -V' * Q * V;
% G = -2 * Q * V;
% 
% function [C,Ceq,G,Geq] = constrfunMarginSAR(Q,V)
% 
% Ceq = [];
% C = reshape(pagemtimes(V, 'transpose', pagemtimes(Q, V), 'none'), size(Q, 3), 1) - 1;
% Geq = [];
% G = 2 * reshape(pagemtimes(Q, V), numel(V), size(Q, 3));
% 
% function val = SAR(Q,V)
% val = pagemtimes(V, 'transpose', pagemtimes(Q, V), 'none');
% val = max(val, [], 'all');
% val = max(val, 0);%for numerical stability
% 




