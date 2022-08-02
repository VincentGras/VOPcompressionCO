function [M, x] = computeMarginSAR(Q, Qvop)

Nc = size(Qvop, 1);
N = size(Q, 3);

if (~isreal(Q) || ~isreal(Qvop))
    
    CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
    [M, xr] = computeMarginSAR(CHtoRS(Q), CHtoRS(Qvop));
    x = xr(1:Nc,:) + 1i * xr(Nc+1:2*Nc, :);
    
else
    
    x = zeros(Nc, N);
    M = zeros(1, N);
    
    
    
    
    parfor i = 1:N
        
        Qi = double(Q(:,:,i));
        Qi = 0.5*(Qi+Qi');
        
        if (size(Qvop, 3)==1)
            
            [eigV, eigD] = eig(Qi-Qvop, 'vector');
            eigD = real(eigD);
            [~, k] = max(eigD);
            x(:,i) = real(eigV(:,k));
            
        else
            
            [eigV, eigD] = eig(Qi, 'vector');
            eigD = real(eigD);
            [~, k] = max(eigD);
            x0 = real(eigV(:,k));
            C = SAR(Qvop, x0)+eps;
            x0 = (0.5 / sqrt(C))*x0;
            
            opt = optimset('Algorithm', 'sqp', 'MaxIter', 1000,  'Display', 'off', 'GradObj', 'on', 'GradConstr', 'on'); 
            
            x(:,i) = fmincon(@(x) objfun(Qi,x), x0, [], [], [], [], [], [], @(x) constrfun(Qvop, x), opt);
            
        end
        
        M(i) = SAR(Qvop,x(:,i))./SAR(Qi,x(:,i)) - 1;
        
        if (mod(i,1000)== 0)
            fprintf('parfor loop : %d / %d; %f\n', i, N, M(i));
        end
        
        
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



