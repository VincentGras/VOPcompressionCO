%% fig 2 : NVOP

figure; 
subplot(1,3,1);
plot_nVOP('ResultsMarch23\HugoNova.mat');
title('a)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
subplot(1,3,2);
plot_nVOP('ResultsMarch23\EllaRapid.mat');
title('b)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
subplot(1,3,3);
plot_nVOP('ResultsMarch23\TheloniusAvanti216tx.mat');
title('c)')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
p = get(gcf, 'Position');
set(gcf, 'Position', [48.2000  291.4000  740.0000  211.2000]);

%%

export_fig ('-transparent', 'figures\Fig1_NVOP.jpg')



%% fig 6 : Rmap and vop locations
load ResultsMarch23\TheloniusAvanti216tx.mat epsil
[~,index] = min(abs(epsil-0.05))
if false
    script2303_VOPCheck('ResultsMarch23\TheloniusAvanti216tx.mat', [1, index, 10], 'ico');
end

load ResultsMarch23\TheloniusAvanti216tx.mat SARc_iCO VOPclassif_iCO epsil
load Avanti216TX\Thelonius\Mask.mat Mask

transf = @(x) x;
%transf = @(x) flip(permute(x, [2 1 3 4]), 3);

%

Rmaps = transf(cat(4, vect2map(Mask, SARc_iCO{1}), ...
    vect2map(Mask, SARc_iCO{index}), ...
    vect2map(Mask, SARc_iCO{end})));
voploc = transf(cat(4, ...
    vect2map(Mask, VOPclassif_iCO(1, :), 0), ...
    vect2map(Mask, VOPclassif_iCO(index, :), 0), ...
    vect2map(Mask, VOPclassif_iCO(end, :), 0)));
for i = 1:3
    voploc(:,:,:,i) = 1./bwdist(voploc(:,:,:,i));
end
voploc = applyMask(transf(Mask), voploc, NaN);

%
figure;
tf = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');

i=2;
k=0;

for j = 1:3
    %for i = 1:3
    k=k+1;
    subplot(3,3,k)
    VARDAS_view_plane(i, max(Rmaps(:,:,:,j), [], i), 'cscale', [0,1]);
    %title(['a'+k-1,')'])
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    box off
    axis off
    colormap('hot');
    if (j==3)
        colorbar;
    end
    %end
end

% vop in sagittal ori
i = 2
for j = 1:3
    
    %for i = 1:3
    k=k+1;
    subplot(3,3,k)
    VARDAS_view_plane(i, max(voploc(:,:,:,j), [], i))
    %title(['a'+k-1,')'])
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    box off
    axis off
    colormap('hot');
    
end

% vop in axial ori
i = 3;
for j = 1:3
    
    %for i = 1:3
    k=k+1;
    subplot(3,3,k)
    VARDAS_view_plane(i, max(voploc(:,:,:,j), [], i), 'cscale', [0, 1])
    %title(['a'+k-1,')'])
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    box off
    axis off
    colormap('hot');
        if (j==3)
        colorbar;
    end
end
%%
export_fig ('-transparent', 'ResultsMarch23\Fig6_Rmap.jpg')


