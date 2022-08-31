load testQ.mat Q10g

Qm = 0.004 * eye(8);

%% test computeVOP_CO


[c1, Qvop1] = computeVOP(@computeCriterionSED, Q10g, Qm);
R1sed = computeCriterionSED(Q10g, Qvop1);
R1sar = computeCriterionSAR(Q10g, Qvop1);

max(R1sed)
max(R1sar)

[c2, Qvop2] = computeVOP(@computeCriterionSAR, Q10g, Qm);
R2sed = computeCriterionSED(Q10g, Qvop2);
R2sar = computeCriterionSAR(Q10g, Qvop2);
