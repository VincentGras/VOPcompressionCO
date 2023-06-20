function f = script_VOPCheckOverestimation(f, subset, algo)

% Get SAR overestimation

load (f, 'epsil', 'Q10g')

if (nargin < 2 || isempty(subset))
    subset = 1:numel(epsil);
end

if (nargin < 3)
    algo = '*';
end

myCluster=parcluster('local');
myCluster.NumWorkers=min(32, myCluster.NumWorkers);
PoolObj = parpool(myCluster,myCluster.NumWorkers);

try
    X = gen_RFshim(Q10g);
    
    if ismember(lower(algo), {'*', 'clu'})
        
        load (f, 'epsil', 'overestimation_CLU',  'VOP_CLU')
        
        if ~exist('overestimation_CLU', 'var')
            overestimation_CLU = cell(size(epsil));
        end
        
        for i = subset
            overestimation_CLU{i} = get_overestimation_statistics(Q10g/50, VOP_CLU{i}, X);
            save (f, '-append', 'overestimation_CLU');
        end
    end
    
    if ismember(lower(algo), {'*', 'icc'})
        
        load (f, 'epsil', 'overestimation_iCC',  'VOP_iCC')
        
        if ~exist('overestimation_iCC', 'var')
            overestimation_iCC = cell(size(epsil));
        end
        
        for i = subset
            overestimation_iCC{i} = get_overestimation_statistics(Q10g, VOP_iCC{i}, X);
            save (f, '-append', 'overestimation_iCC');
        end
    end
    
    if ismember(lower(algo), {'*', 'ico'})
        
        load (f, 'epsil', 'overestimation_iCO',  'VOP_iCO')
        
        if ~exist('overestimation_iCO', 'var')
            overestimation_iCO = cell(size(epsil));
        end
        
        for i = subset
            overestimation_iCO{i} = get_overestimation_statistics(Q10g, VOP_iCO{i}, X);
            save (f, '-append', 'overestimation_iCO');
        end
    end
end

delete(PoolObj)

function X = gen_RFshim(Q10g)
nX = 1e5;
Nc = size(Q10g, 1);
X = randn(Nc,nX) + 1i * randn(Nc,nX);
X = X ./ vecnorm(X, 2, 1);

function [stat] = get_overestimation_statistics(Q, Qvop, X)

Nc = size(X, 1);
nX = size(X, 2);
assert(Nc == size(Q, 1), 'nc = %d ???', Nc);
overest = zeros(nX, 1);

parfor i = 1:nX
    SAR_Q_ = pagemtimes(X(:,i), 'ctranspose', pagemtimes(Q, 'none', X(:,i), 'none'), 'none');
    SAR_Qvop_ = pagemtimes(X(:,i), 'ctranspose', pagemtimes(Qvop, 'none', X(:,i), 'none'), 'none');
    overest(i) = max(real(SAR_Qvop_)) / max(real(SAR_Q_));
end

stat = [min(overest), median(overest), max(overest), mean(overest)];






