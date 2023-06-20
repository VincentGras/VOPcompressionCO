function f = script_VOPCheck(f, subset, algo)

% run rQstar(VOP, Q) on each SAR matrices Q of the SAR model save in file f

load (f, 'epsil', 'Q10g')

if (nargin < 2 || isempty(subset))
    subset = 1:numel(epsil);
end

if (nargin < 3)
    algo = '*';
end

criterion = @(Q,Qvop) rQstar(Q, Qvop);

myCluster=parcluster('local');
myCluster.NumWorkers=min(48, myCluster.NumWorkers);
PoolObj = parpool(myCluster, myCluster.NumWorkers);

try
    if ismember(lower(algo), {'*', 'clu'})
        
        load (f, 'epsil', 'SARc_CLU',  'VOP_CLU')
        
        if ~exist('SARc_CLU', 'var')
            SARc_CLU = cell(size(epsil));
        end
        if ~exist('ct_vopcheck_CLU', 'var')
            ct_vopcheck_CLU = zeros(size(epsil));
        end
        
        for i = subset
            tic;
            SARc_CLU{i} = feval(criterion, Q10g/50, VOP_CLU{i});
            ct_vopcheck_CLU(i) = toc;
            save (f, '-append', 'SARc_CLU', 'ct_vopcheck_CLU');
        end
    end
    
    if ismember(lower(algo), {'*', 'icc'})
        
        load (f, 'epsil', 'SARc_iCC',  'VOP_iCC')
        
        if ~exist('SARc_iCC', 'var')
            SARc_iCC = cell(size(epsil));
        end
        if ~exist('ct_vopcheck_iCC', 'var')
            ct_vopcheck_iCC = zeros(size(epsil));
        end
        
        for i = subset
            tic;
            SARc_iCC{i} = feval(criterion, Q10g, VOP_iCC{i});
            ct_vopcheck_iCC(i) = toc;
            save (f, '-append', 'SARc_iCC', 'ct_vopcheck_iCC');
        end
    end
    
    if ismember(lower(algo), {'*', 'ico'})
        
        load (f, 'epsil', 'SARc_iCO',  'VOP_iCO')
        
        
        if ~exist('SARc_iCO', 'var')
            SARc_iCO = cell(size(epsil));
        end
        
        if ~exist('ct_vopcheck_iCO', 'var')
            ct_vopcheck_iCO = zeros(size(epsil));
        end
        
        for i = subset
            tic;
            SARc_iCO{i} = feval(criterion, Q10g, VOP_iCO{i});
            ct_vopcheck_iCO(i) = toc;
            save (f, '-append', 'SARc_iCO', 'ct_vopcheck_iCO');
        end
    end
end

delete(PoolObj)





