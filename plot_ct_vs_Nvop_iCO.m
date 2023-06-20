function plot_ct_vs_Nvop_iCO(simu, newfig)

if (nargin < 2)
    newfig=false;
end

load (simu, 'ct_iCO', 'nVOP_iCO')

L = {'iCO full model'};

if newfig
    figure;
end

loglog(    nVOP_iCO, cumsum(ct_iCO)/60,'o--', 'Color', 'b', 'MarkerFaceColor', 'b');

xlabel('N^*')
%xlim([0, 25]);
ylabel('Computation time (min)');
box on;
grid on;
legend(L);




