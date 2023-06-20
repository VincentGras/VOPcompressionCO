%% Fig 4 : CT against N* (log-log plot)

fprintf('TheloniusAvanti216tx\n')
load ResultsMarch23\TheloniusAvanti216tx.mat nVOP_CLU ct_CLU nVOP_iCC ct_iCC nVOP_iCO ct_iCO
data_Afull = ct_logloganalysis(nVOP_CLU, ct_CLU, nVOP_iCC, ct_iCC, nVOP_iCO, ct_iCO);

fprintf('TheloniusAvanti216tx_subsamp8\n')
load ResultsMarch23\TheloniusAvanti216tx_subsamp8.mat nVOP_CLU ct_CLU nVOP_iCC ct_iCC nVOP_iCO ct_iCO
data_A = ct_logloganalysis(nVOP_CLU, ct_CLU, nVOP_iCC, ct_iCC, nVOP_iCO, ct_iCO);

fprintf('HugoNova\n')
load ResultsMarch23\HugoNova.mat nVOP_CLU ct_CLU nVOP_iCC ct_iCC nVOP_iCO ct_iCO
data_N = ct_logloganalysis(nVOP_CLU, ct_CLU, nVOP_iCC, ct_iCC, nVOP_iCO, ct_iCO);

fprintf('EllaRapid\n')
load ResultsMarch23\EllaRapid.mat nVOP_CLU ct_CLU nVOP_iCC ct_iCC nVOP_iCO ct_iCO
data_R = ct_logloganalysis(nVOP_CLU, ct_CLU, nVOP_iCC, ct_iCC, nVOP_iCO, ct_iCO);
%%

figure;
subplot(1,3,1);
loglog(data_N.nVOP_CLU, data_N.ct_CLU/60, 'ko', ...
    data_N.nVOP_iCC, data_N.ct_iCC/60, 'rx', ...
    data_N.nVOP_iCO, data_N.ct_iCO/60, 'bs', ...
    data_N.nVOP_CLU, exp(polyval(data_N.ab_CLU, log(data_N.nVOP_CLU)))/60, 'k-', ...
    data_N.nVOP_iCC, exp(polyval(data_N.ab_iCC, log(data_N.nVOP_iCC)))/60, 'r-', ...
    data_N.nVOP_iCO, exp(polyval(data_N.ab_iCO, log(data_N.nVOP_iCO)))/60, 'b-');
xlabel('N^*')
ylabel('CT (min)')
ylim([0.005,4e2])

legend(sprintf('CLU (%.1f)', data_N.ab_CLU(1)), ...
       sprintf('iCC (%.1f)', data_N.ab_iCC(1)), ...
       sprintf('iCO (%.1f)', data_N.ab_iCO(1)), 'Location', 'SouthEast');
grid on
title('a)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';


subplot(1,3,2);
loglog(data_R.nVOP_CLU, data_R.ct_CLU/60, 'ko', ...
    data_R.nVOP_iCC, data_R.ct_iCC/60, 'rx', ...
    data_R.nVOP_iCO, data_R.ct_iCO/60, 'bs', ...
    data_R.nVOP_CLU, exp(polyval(data_R.ab_CLU, log(data_R.nVOP_CLU)))/60, 'k-', ...
    data_R.nVOP_iCC, exp(polyval(data_R.ab_iCC, log(data_R.nVOP_iCC)))/60, 'r-', ...
    data_R.nVOP_iCO, exp(polyval(data_R.ab_iCO, log(data_R.nVOP_iCO)))/60, 'b-');
ylim([0.005,1e2])
xlabel('N^*')
ylabel('CT (min)')
legend(sprintf('CLU (%.1f)', data_R.ab_CLU(1)), ...
       sprintf('iCC (%.1f)', data_R.ab_iCC(1)), ...
       sprintf('iCO (%.1f)', data_R.ab_iCO(1)), 'Location', 'SouthEast');
grid on
title('b)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

subplot(1,3,3);
loglog(data_Afull.nVOP_CLU, data_Afull.ct_CLU/60, 'ko', ...
    data_Afull.nVOP_iCC, data_Afull.ct_iCC/60, 'rx', ...
    data_Afull.nVOP_iCO, data_Afull.ct_iCO/60, 'bs', ...
    data_Afull.nVOP_CLU, exp(polyval(data_Afull.ab_CLU, log(data_Afull.nVOP_CLU)))/60, 'k-', ...
    data_Afull.nVOP_iCC, exp(polyval(data_Afull.ab_iCC, log(data_Afull.nVOP_iCC)))/60, 'r-', ...
    data_Afull.nVOP_iCO, exp(polyval(data_Afull.ab_iCO, log(data_Afull.nVOP_iCO)))/60, 'b-');
xlabel('N^*')
ylabel('CT (min)')
ylim([0.01,5e4])
legend(sprintf('CLU (%.1f)', data_Afull.ab_CLU(1)), ...
       sprintf('iCC (%.1f)', data_Afull.ab_iCC(1)), ...
       sprintf('iCO (%.1f)', data_Afull.ab_iCO(1)), 'Location', 'SouthEast');
grid on
title('c)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

%%

p = get(gcf, 'Position');
set(gcf, 'Position', [12.2000   69.0000  740.0000  211.2000]);
export_fig ('-transparent', 'Fig3_CT_loglog.jpg')

