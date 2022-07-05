function S = spectralNorm(Q)

N = size(Q, 3); 

S = zeros(N,1);

parfor i = 1:N
    S(i) = norm(Q(:,:,i), 2);
end
