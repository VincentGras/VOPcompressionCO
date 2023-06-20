%% Fig 6 : CT against N* (log-log plot)
% !! run script_VOPcheck first !!

fprintf('TheloniusAvanti216tx\n')
load ResultsMarch23\TheloniusAvanti216tx.mat nVOP_CLU ct_vopcheck_CLU nVOP_iCC ct_vopcheck_iCC nVOP_iCO ct_vopcheck_iCO
ab_Afull = ct_vopcheck_logloganalysis(nVOP_iCO, ct_vopcheck_iCO);
nVOP_iCO_Afull = nVOP_iCO;
ct_vopcheck_iCO_Afull = ct_vopcheck_iCO;

fprintf('TheloniusAvanti216tx_subsamp8\n')
load ResultsMarch23\TheloniusAvanti216tx_subsamp8.mat nVOP_CLU ct_vopcheck_CLU nVOP_iCC ct_vopcheck_iCC nVOP_iCO ct_vopcheck_iCO
ab_A = ct_vopcheck_logloganalysis(nVOP_iCO, ct_vopcheck_iCO);
nVOP_iCO_A = nVOP_iCO;
ct_vopcheck_iCO_A = ct_vopcheck_iCO;

fprintf('HugoNova\n')
load ResultsMarch23\HugoNova.mat nVOP_CLU ct_vopcheck_CLU nVOP_iCC ct_vopcheck_iCC nVOP_iCO ct_vopcheck_iCO
ab_N = ct_vopcheck_logloganalysis(nVOP_iCO, ct_vopcheck_iCO);
nVOP_iCO_N = nVOP_iCO;
ct_vopcheck_iCO_N = ct_vopcheck_iCO;

fprintf('EllaRapid\n')
load ResultsMarch23\EllaRapid.mat nVOP_CLU ct_vopcheck_CLU nVOP_iCC ct_vopcheck_iCC nVOP_iCO ct_vopcheck_iCO
ab_R = ct_vopcheck_logloganalysis(nVOP_iCO, ct_vopcheck_iCO);
nVOP_iCO_R = nVOP_iCO;
ct_vopcheck_iCO_R = ct_vopcheck_iCO;

%%

figure;
loglog(nVOP_iCO_N, ct_vopcheck_iCO_N/60, 'b*', ... 
    nVOP_iCO_R, ct_vopcheck_iCO_R/60, 'ro', ... 
    nVOP_iCO_Afull, ct_vopcheck_iCO_Afull/60, 'kx', ...
    nVOP_iCO_N, exp(polyval(ab_N, log(nVOP_iCO_N)))/60, 'b-', ...
    nVOP_iCO_R, exp(polyval(ab_R, log(nVOP_iCO_R)))/60, 'r-', ...
    nVOP_iCO_Afull, exp(polyval(ab_Afull, log(nVOP_iCO_Afull)))/60, 'k-');
xlabel('N^*')
ylabel('CT(r^*_{max}) (min)')
legend(sprintf('Nova (%.1f)', ab_N(1)), ...
       sprintf('Rapid (%.1f)', ab_R(1)), ...
       sprintf('Avanti (%.1f)', ab_Afull(1)), 'Location', 'NorthWest');
grid on;

%%

p = get(gcf, 'Position');
set(gcf, 'Position', [12.2000   69.0000  400 350]);
export_fig ('-transparent', 'Fig5_CTcheck_loglog.jpg')
