%% Fig 6 : SAR overestimation (\eta) versus N*
% !! run script_VOPcheckOverestimation first !!

fprintf('TheloniusAvanti216tx_subsamp8\n')
load ResultsMarch23\TheloniusAvanti216tx.mat nVOP_CLU overestimation_CLU nVOP_iCC overestimation_iCC nVOP_iCO overestimation_iCO epsil

nVOP_CLU_A = nVOP_CLU;
overestimation_CLU_A = cat(1, overestimation_CLU{:});
nVOP_iCC_A = nVOP_iCC;
overestimation_iCC_A = cat(1, overestimation_iCC{:});
nVOP_iCO_A = nVOP_iCO;
overestimation_iCO_A = cat(1, overestimation_iCO{:});

fprintf('HugoNova\n')
load ResultsMarch23\HugoNova.mat nVOP_CLU overestimation_CLU nVOP_iCC overestimation_iCC nVOP_iCO overestimation_iCO
nVOP_CLU_N = nVOP_CLU;
overestimation_CLU_N = cat(1, overestimation_CLU{:});
nVOP_iCC_N = nVOP_iCC;
overestimation_iCC_N = cat(1, overestimation_iCC{:});
nVOP_iCO_N = nVOP_iCO;
overestimation_iCO_N = cat(1, overestimation_iCO{:});

fprintf('EllaRapid\n')
load ResultsMarch23\EllaRapid.mat nVOP_CLU overestimation_CLU nVOP_iCC overestimation_iCC nVOP_iCO overestimation_iCO
nVOP_CLU_R = nVOP_CLU;
overestimation_CLU_R = cat(1, overestimation_CLU{:});
nVOP_iCC_R = nVOP_iCC;
overestimation_iCC_R = cat(1, overestimation_iCC{:});
nVOP_iCO_R = nVOP_iCO;
overestimation_iCO_R = cat(1, overestimation_iCO{:});
%%
figure; 
subplot(1,3,1);
semilogx(nVOP_CLU_N(1:size(overestimation_CLU_N,1)), overestimation_CLU_N(:, end), 'ko-', ...
    nVOP_iCC_N, overestimation_iCC_N(:, end), 'rx-', ...
    nVOP_iCO_N, overestimation_iCO_N(:, end), 'bs-')
xlabel('N^*')
ylabel('\eta')
legend('CLU', 'iCC', 'iCO', 'Location', 'NorthEast');
grid on
title('a)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
subplot(1,3,2);
semilogx(nVOP_CLU_R(1:size(overestimation_CLU_R,1)), overestimation_CLU_R(:, end), 'ko-', ...
    nVOP_iCC_R, overestimation_iCC_R(:, end), 'rx-', ...
    nVOP_iCO_R, overestimation_iCO_R(:, end), 'bs-')
xlabel('N^*')
ylabel('\eta')
legend('CLU', 'iCC', 'iCO', 'Location', 'NorthEast');
grid on
title('b)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
subplot(1,3,3);
semilogx(nVOP_CLU_A(1:size(overestimation_CLU_A,1)), overestimation_CLU_A(:, end), 'ko-', ...
    nVOP_iCC_A(1:size(overestimation_iCC_A,1)), overestimation_iCC_A(:, end), 'rx-', ...
    nVOP_iCO_A(1:size(overestimation_iCO_A,1)), overestimation_iCO_A(:, end), 'bs-')
xlabel('N^*')
ylabel('\eta')
legend('CLU', 'iCC', 'iCO', 'Location', 'NorthEast');
grid on
title('c)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

p = get(gcf, 'Position');
set(gcf, 'Position', [12.2000   69.0000  740.0000  211.2000]);

%%

export_fig ('-transparent', 'ResultsMarch23\Fig5_overestimation_loglog.jpg')
