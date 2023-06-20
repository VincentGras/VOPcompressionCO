%% fig 3 : CT

figure; 
subplot(1,3,1);
plot_ct('ResultsMarch23\HugoNova.mat');
title('a)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
subplot(1,3,2);
plot_ct('ResultsMarch23\EllaRapid.mat');
title('b)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
subplot(1,3,3);
plot_ct('ResultsMarch23\TheloniusAvanti216tx.mat');
title('b)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
p = get(gcf, 'Position');
set(gcf, 'Position', [12.2000   69.0000  740.0000  211.2000]);

%%

export_fig ('-transparent', 'figures\Fig2_CT.jpg')


