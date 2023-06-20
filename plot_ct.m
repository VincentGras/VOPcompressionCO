function plot_ct(simu, newfig)

if (nargin < 2)
    newfig=false;
end

load (simu, 'epsil', 'ct_CLU', 'ct_iCC', 'ct_iCO', 'nVOP_CLU', 'nVOP_iCC', 'nVOP_iCO')

ct_iCC(nVOP_iCC == 301) = NaN;

L = {'CLU', 'iCC', 'iCO'};

if newfig
    figure;
end

loglog(epsil * 100, ct_CLU/60, 'o-', 'Color', 'k');
hold on;
semilogy(    epsil * 100, cumsum(ct_iCC)/60,'x-', 'Color', 'r');
semilogy(    epsil * 100, cumsum(ct_iCO)/60,'s--', 'Color', 'b');

xlabel('\lambda'' (%)')
xlim([0, 25]);
ylabel('CT (min)');
box on;
grid on;
legend(L);