function [CSAR, x] = computeCriterionSAR2(Q, Qvop, cutoff, r)


Nc = size(Qvop, 1);
N = size(Q, 3);

if (nargin < 4 || isinf(r) || isnan(r))
    r = computeCriterionSAR(eye(Nc), Qvop);
    fprintf('********* r = %f\n', r)
end

if (~isreal(Q) || ~isreal(Qvop))
    
    CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
    [CSAR, xr] = computeCriterionSAR2(CHtoRS(Q), CHtoRS(Qvop), cutoff, r);
    x = xr(1:Nc,:) + 1i * xr(Nc+1:2*Nc, :);
    
else
    
    x = zeros(Nc, N);
    CSAR = zeros(1, N);
    
    parfor i = 1:N
        
        Qi = double(Q(:,:,i));
        Qi = 0.5*(Qi+Qi');
        
        if (size(Qvop, 3)==1)
            
            [eigV, eigD] = eig(Qi-Qvop, 'vector');
            eigD = real(eigD);
            [~, k] = max(eigD);
            x(:,i) = real(eigV(:,k));
            CSAR(i) = SAR(Qi,x(:,i))./SAR(Qvop,x(:,i));
        else
            
           

            % SVD analysis 
            [U, eigD] = eig(Qi, 'vector');
            eigD = real(eigD);
            [eigD, k] = sort(abs(eigD), 'descend');
            U = real(U(:,k));

            E = eigD>cutoff/r; %*eigD(1);
            Edim = nnz(E);
            
            if (Edim > 0)
                %F = ~E;
                %Fdim = nnz(F);
                QiSVD = diag(eigD(E));

                QvopSVD = pagemtimes(Qvop, U);
                QvopSVD = pagemtimes(U', QvopSVD);

                QvopSVD_E = QvopSVD(E, E, :); % m*m
                %QvopSVD_EF = QvopSVD(E, F, :); % m*n

                x0 = zeros(Edim, 1); 
                x0(1) = 1;
                x0 = sqrt(0.5 / (SAR(QvopSVD_E, x0)+eps))*x0;

                opt = optimset('Algorithm', 'sqp', 'MaxIter', 1000,  'Display', 'off', 'GradObj', 'on', 'GradConstr', 'on'); 

                
                yiE = fmincon(@(x) objfun(QiSVD,x), x0, [], [], [], [], [], [], @(x) constrfun(QvopSVD_E, x), opt);
                yi = zeros(Nc,1);
                yi(E) = yiE;
                xi = U * yi;
                x(:, i) = xi;
                
                %[ SAR(QiSVD,yiE)./SAR(QvopSVD_E,yiE), SAR(Qi,x(:,i))./SAR(Qvop,x(:,i))]
               
                CSAR(i) = cutoff  + SAR(QiSVD,yiE)./SAR(QvopSVD_E,yiE);
            end
        end
        
        
        
        if (mod(i,1000)== 0)
            fprintf('parfor loop : %d / %d; %f\n', i, N, CSAR(i));
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



