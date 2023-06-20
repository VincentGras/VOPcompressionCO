%% fig 6 : rQstar maps and vop locations in MIP

load ResultsMarch23\TheloniusAvanti216tx.mat SARc_iCO VOPclassif_iCO epsil
load Avanti216TX\Thelonius\Mask.mat Mask
[~,index] = min(abs(epsil - 0.05));

Rmaps = cat(4, vect2map(Mask, SARc_iCO{1}), ...
    vect2map(Mask, SARc_iCO{index}), ...
    vect2map(Mask, SARc_iCO{end}));
voploc = cat(4, ...
    vect2map(Mask, VOPclassif_iCO(1, :), 0), ...
    vect2map(Mask, VOPclassif_iCO(index, :), 0), ...
    vect2map(Mask, VOPclassif_iCO(end, :), 0));
for i = 1:3
    voploc(:,:,:,i) = 1./bwdist(voploc(:,:,:,i));
end
voploc = applyMask(Mask, voploc, NaN);

%% r_Q* maps

figure;
tf = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');

i = 2; % sagittal view
k = 0;

for j = 1:3 % loop across 3 different values of compression parameters
    k=k+1;
    subplot(3,3,k)
    VARDAS_view_plane(i, max(Rmaps(:,:,:,j), [], i), 'cscale', [0,1]);
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    box off
    axis off
    colormap('hot');
    if (j==3)
        colorbar;
    end
end

% vop localizations in sagittal ori
i = 2; 
for j = 1:3
    k=k+1;
    subplot(3,3,k)
    VARDAS_view_plane(i, max(voploc(:,:,:,j), [], i))
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    box off
    axis off
    colormap('hot');    
end

% vop localizations in axial ori
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

%% expo

export_fig ('-transparent', 'Fig6_Rmap.jpg')

