function plot_ct_vs_Nvop(simu, newfig)

if (nargin < 2)
    newfig=false;
end

load (simu, 'ct_CLU', 'ct_iCC', 'ct_iCO', 'nVOP_CLU', 'nVOP_iCC', 'nVOP_iCO')

ct_iCC(nVOP_iCC == 301) = NaN;

L = {'CLU', 'iCC', 'iCO'};

if newfig
    figure;
end

loglog(nVOP_CLU, ct_CLU/60, 'o-', 'Color', 'k');
hold on;
loglog(    nVOP_iCC, cumsum(ct_iCC)/60,'x-', 'Color', 'r');
loglog(    nVOP_iCO, cumsum(ct_iCO)/60,'s--', 'Color', 'b');

xlabel('N^*')
%xlim([0, 25]);
ylabel('Computation time (min)');
box on;
grid on;
legend(L);




