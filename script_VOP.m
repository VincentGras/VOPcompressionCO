function outputfile = script2303_VOP(choice, testmyscript, subsamp )

% Launch VOP compression
% choice = nova | rapidElla | avanti216tx
% testmyscript = test | [notAtest]
% subsamp = > 0 integer [1]

addpath CO
addpath CC
addpath CLU\

if (nargin == 0)
    myCluster=parcluster('local');
    myCluster.NumWorkers=min(48, myCluster.NumWorkers);
    PoolObj = parpool(myCluster,myCluster.NumWorkers);
    try
        script2303_VOP nova notAtest
    end
    try
        script2303_VOP rapidella notAtest
    end
    try
        script2303_VOP avanti216tx notAtest
    end
    delete(PoolObj)
    return;
end

if (nargin < 3)
    subsamp = 1;
end

if (ischar(subsamp))
    subsamp = str2num(subsamp);
end

if (nargin > 1 && ischar(testmyscript) && strcmp(testmyscript,'test'))
    testmyscript = true;
else
    testmyscript = false;
    
end

Nc=8;
switch (lower(choice))
    case 'nova'
        inputfile_Mask = 'Nova/Hugo/Mask_HeadModel.mat';
        inputfile_Q = 'Nova/Hugo/Q10gMatrices_HeadModel.mat';
        outputfile = 'HugoNova.mat';
    case 'rapid'
        inputfile_Mask = 'Rapid/Giorgio/Mask_Giorgio.mat';
        inputfile_Q = 'Rapid/Giorgio/Q10gMatrices_Giorgio.mat';
        outputfile = 'GiorgioRapid.mat';
    case 'rapidella'
        inputfile_Mask = 'Rapid/Ella/Mask_Ella.mat';
        inputfile_Q = 'Rapid/Ella/Q10gMatrices_Ella.mat';
        outputfile = 'EllaRapid.mat';
    case 'avanti216tx'
        inputfile_Mask = 'Avanti216TX/Thelonius/Mask.mat';
        inputfile_Q = 'Avanti216TX/Thelonius/Q_16TX.mat';
        outputfile = 'TheloniusAvanti216tx.mat';
        
        if (subsamp > 1)
            outputfile = sprintf('TheloniusAvanti216tx_subsamp%d.mat', subsamp);
        end

        Nc = 16;
    otherwise
        error ('cas inconnu');
end


% if testmyscript
%     subsamp = subsamp * 1000;
% end

if (testmyscript)
    outputfile = ['test__', outputfile];
end

if (~isdir('ResultsMarch23'))
    mkdir('ResultsMarch23')
end

outputfile = fullfile('ResultsMarch23', outputfile);


%%

%load(inputfile_Mask, 'Mask')
load(inputfile_Q, 'Q10g', 'QGlobal')
if (size(Q10g, 1) ~= Nc)
    Q10g = permute(Q10g, [2 3 1]);
end
Q10g=double(Q10g*50);
QGlobal = double(QGlobal*50);

%% susbampling and enforce Hermitian symmetry

Q10g = Q10g(:,:,1:subsamp:end);
Q10g = 0.5 * (Q10g + conj( permute(Q10g, [2 1 3])));

%% max(|Q10g|) ?

normQ10g = spectralNorm(Q10g);
normQ10gmax = max(normQ10g);
normQGlobal = norm(QGlobal);

save (outputfile, 'Q10g', 'QGlobal', 'normQ10gmax', 'normQGlobal');


%% SAR model compression settings

maxnvop = 400;
iterparam = 0.75;
epsilmin = 0.015;
niter = 10;
epsilmax =  0.9999 * epsilmin * (1/iterparam)^(niter-1);
epsil = epsilmax * iterparam.^(0:niter-1);
save (outputfile, '-append', 'epsil', 'epsilmax', 'niter', 'epsilmin', 'iterparam', 'maxnvop');

%% iCO

epsil_iCO  = epsil;
[VOPclassif_iCO, VOP_iCO, ~, ct_iCO] = ...
    computeVOPi_CO(Q10g, epsilmax * normQ10gmax, niter, iterparam, maxnvop, []);
nVOP_iCO = sum(VOPclassif_iCO==1, 2);
save (outputfile, '-append', 'epsil_iCO', 'VOP_iCO', 'nVOP_iCO', 'ct_iCO', 'VOPclassif_iCO');

figure;
plot(epsil, nVOP_iCO)
drawnow

criterion = @(Q,Qvop) rQstar(Q, Qvop);

SARc_iCO = verif(criterion, Q10g, VOP_iCO);
save (outputfile, '-append', 'SARc_iCO');

%% iCCC

launchtime = datestr(now);
launchtime = regexprep(launchtime, '\s', '_');
launchtime = regexprep(launchtime, ':', '');

if (~isdir('temp'))
    mkdir('temp');
end

fname_in_VOP_Compression = [choice, '_', launchtime, '.mat'];
fname_in_VOP_Compression = fullfile('temp', fname_in_VOP_Compression);

matrices = Q10g;
save (fname_in_VOP_Compression,  'matrices');
clear matrices

[oiCCC] = ...
    VOP_compression_iterative_VG(fname_in_VOP_Compression, epsilmax, iterparam, eye(Nc), choice, 0, maxnvop, epsilmin);

niter_iCC = numel(oiCCC);

epsil_iCC = zeros(niter_iCC, 1);
VOP_iCC = cell(niter_iCC, 1);
nVOP_iCC = zeros(niter_iCC, 1);
ct_iCC = zeros(niter_iCC, 1);

for i = 1:niter_iCC
    data = load(oiCCC{i});
    epsil_iCC(i) = data.eps_G;
    VOP_iCC{i} = data.VOP;
    nVOP_iCC(i) = size(data.VOP, 3);
    ct_iCC(i) = data.elapsed_time;
end

save (outputfile, '-append', 'iterparam',  'maxnvop', 'niter_iCC', 'epsil_iCC', 'VOP_iCC', 'nVOP_iCC', 'ct_iCC');

hold on;
plot(epsil_iCC, nVOP_iCC);
drawnow;

SARc_iCC = verif(criterion, Q10g, VOP_iCC);
save (outputfile, '-append', 'SARc_iCC');

%% CLU

epsil_CLU = epsil;
save (outputfile, '-append', 'epsil_CLU');
VOP_CLU = cell(niter, 1);
ct_CLU = NaN(niter, 1);
nVOP_CLU = NaN(niter, 1);


for i = 1:niter

    % CLU
    
    if (i == 1 || nVOP_CLU(i-1) < maxnvop)
        try
            tic;
            VOP_CLU{i} = VOPcompressionCLU(fullfile('VOPcompressionCLU.dat'), ...
                Q10g/50, epsil(i));
            ct_CLU(i) = toc;
            nVOP_CLU(i) = size(VOP_CLU{i}, 3);
            save (outputfile, '-append', 'VOP_CLU', 'nVOP_CLU', 'ct_CLU');
        end
    end
    
end

hold on;
plot(epsil_CLU, nVOP_CLU)
title('NVOP')
xlabel('epsil')
legend('iCO', 'iCC', 'CLU');

figure; 
plot(epsil_iCO, ct_iCO, epsil_iCC, ct_iCC, epsil_CLU, ct_CLU);
title('CT (s)')
xlabel('epsil')
legend('iCO', 'iCC', 'CLU');

SARc_CLU = verif(criterion, Q10g/50, VOP_CLU);
save (outputfile, '-append', 'SARc_CLU');

function v = verif(metric, Q10g, Qvop)
v = cell(size(Qvop));
for i = 1:numel(v)
    try
        v{i} = feval(metric, Q10g, Qvop{i});
    end
end




