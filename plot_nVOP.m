function plot_nVOP(simu)

load (simu, 'epsil', 'nVOP_CLU', 'nVOP_iCC', 'nVOP_iCO')

nVOP_iCC(nVOP_iCC == 301) = NaN;

L = {'CLU', 'iCC', 'iCO'};

loglog( epsil * 100, (nVOP_CLU), 'o-', 'Color', 'k'); %, 'MarkerFaceColor', 'b');
yticks([1, 50, 100, 300, 800])
hold on; 
loglog( epsil * 100, (nVOP_iCC), 'x-', 'Color', 'r');
loglog( epsil * 100, (nVOP_iCO), 's--', 'Color', 'b');

    xlabel('\lambda'' (%)')
xlim([0, 25]);
%ylim([5,500]);
ylabel('N^*');
box on;
grid on;
legend(L);
