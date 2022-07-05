function [R, V] = testQmatrixDomination_CO(Q, Qvop)

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

N = size(Q, 3);
Nc = size(Qvop, 1);

V = zeros(Nc, N);
R = zeros(1, N);

CHtoRS = @(Q) cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
QvopR = CHtoRS(Qvop);

evalSARReal = @(Q,V) max(max(reshape(pagemtimes(V, 'transpose', pagemtimes(Q, V), 'none'), size(Q, 3), 1)),0);

parfor i = 1:N
    
    QtR = CHtoRS(Q(:,:,i));
    QtR = 0.5*(QtR+QtR');
    
    if (size(Qvop, 3)==1)
        
        [eigV, eigD] = eig(QtR-QvopR, 'vector');
        eigD = real(eigD);
        [~, k] = max(eigD);
        VR = real(eigV(:,k));
        
    else
        
        [eigV, eigD] = eig(QtR, 'vector');
        eigD = real(eigD);
        [~, k] = max(eigD);
        V0R = real(eigV(:,k));
        C = max(evalSARReal(QvopR, V0R));
        V0R = 0.5*V0R / (sqrt(C)+eps);
        
        opt = optimset('Algorithm', 'sqp', 'MaxIter', 1000,  'Display', 'off', 'GradObj', 'on', 'GradConstr', 'on'); %,'OutputFcn',outfun);
        
        VR = fmincon(@(VR) objfun(QtR,VR), V0R, [], [], [], [], [], [], @(VR) VopConstr(QvopR, VR), opt);
    end
    
    V(:,i) = VR(1:Nc) + 1i * VR(Nc+1:2*Nc);
    R(i) = evalSARReal(QtR,VR)./evalSARReal(QvopR,VR);
    
    if (mod(i,1000)== 0)
        fprintf('parfor loop : %d / %d; %f\n', i, N, R(i));
    end
    
    
end

function [C, G] = objfun(Q, V)
C = -V' * Q * V;
G = -2 * Q * V;


function [C,Ceq,G,Geq] = VopConstr(Q,V)

Ceq = [];
C = reshape(pagemtimes(V, 'transpose', pagemtimes(Q, V), 'none'), size(Q, 3), 1) - 1;
Geq = [];
G = 2 * reshape(pagemtimes(Q, V), numel(V), size(Q, 3));



