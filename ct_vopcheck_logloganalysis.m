function ab = ct_vopcheck_logloganalysis(nVOP, ct_vopcheck)

%ab_iCC = polyfit( log(nVOP_iCC), log(ct_vopcheck_iCC), 1);
ab = polyfit( log(nVOP), log(ct_vopcheck), 1);
fprintf('a=%10.2f b=%10.2f\n', ab(1), ab(2));

% figure; 
% loglog(nVOP, ct_vopcheck, 'b'); 
% hold on;
% loglog(nVOP, exp(polyval(ab, log(nVOP))), 'b--');

