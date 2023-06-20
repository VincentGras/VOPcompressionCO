function data = ct_logloganalysis(nVOP_CLU, ct_CLU, nVOP_iCC, ct_iCC, nVOP_iCO, ct_iCO)

data.ab_CLU = polyfit( log(nVOP_CLU(nVOP_CLU>0)), log(ct_CLU(nVOP_CLU>0)), 1);
data.ab_iCC = polyfit( log(nVOP_iCC), log(ct_iCC), 1);
data.ab_iCO = polyfit( log(nVOP_iCO), log(ct_iCO), 1);
fprintf('CLU : a=%10.2f b=%10.2f\n', data.ab_CLU(1), data.ab_CLU(2));
fprintf('iCC : a=%10.2f b=%10.2f\n', data.ab_iCC(1), data.ab_iCC(2));
fprintf('iCO : a=%10.2f b=%10.2f\n', data.ab_iCO(1), data.ab_iCO(2));

data.nVOP_CLU = nVOP_CLU;
data.ct_CLU = ct_CLU;
data.nVOP_iCC = nVOP_iCC;
data.ct_iCC = ct_iCC;
data.nVOP_iCO = nVOP_iCO;
data.ct_iCO = ct_iCO;


% figure; 
% loglog(nVOP_iCC, ct_iCC, 'k'); 
% hold on;
% loglog(nVOP_iCO, ct_iCO, 'r'); 
% loglog(nVOP_iCC, exp(polyval(ab_iCC, log(nVOP_iCC))), 'k--');
% loglog(nVOP_iCO, exp(polyval(ab_iCO, log(nVOP_iCO))), 'r--');
