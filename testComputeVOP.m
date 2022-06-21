load testQ.mat Q10g

Qm = 0.004 * eye(8);

%% test computeVOP_CO

[Qvop1, c1] = computeVOP_CO(Q10g, Qm);
R1 = testQmatrixDomination_CO(Q10g, Qvop1);
max(R1)

%% test computeVOP_iCO

[Qvop2, c2] = computeVOP_iCO(Q10g, Qm, 5, 0.8);
R2 = testQmatrixDomination_CO(Q10g, Qvop2);
max(R2)

%% plot R1, R2

figure; 
plot(R1);
hold on
plot(R2)

%% Verif with random Xs

N = 1e2;
Nc = size(Q10g, 1);
X = exp (1i * (2*pi*rand(Nc, N)));

SAR = zeros(N,1);
SARvop = zeros(N,1);

parfor i = 1:N
    SAR(i) = max(real(pagemtimes(X(:,i), 'ctranspose', pagemtimes(Q10g,X(:,i)), 'none'))); 
    SARvop(i) = max(real(pagemtimes(X(:,i), 'ctranspose', pagemtimes(Qvop2,X(:,i)), 'none'))); 
end

figure;
plot(SAR, SARvop, 'o');
xlabel('SAR')
ylabel('VOP-SAR')
hold on;
plot([0,1],[0,1], 'r', 'LineWidth', 2);
box on;
grid on;
title ('VOP-SAR versus SAR');

rho = SARvop./SAR;

figure; 
histogram(rho, 'Normalization', 'pdf', 'DisplayStyle', 'stairs');
box on;
grid on;
xlabel('rho')
ylabel('count')
title('SAR overstimation')


