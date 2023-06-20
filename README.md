Here is a decription of the content of this code repository :

CO/computeVOP_CO
    VOP compression algorithm based on convex optimization
CO/computeVOPi_CO
	Iterative application of computeVOP_CO with dicreasing compression margin following the approach described by Orzada er al (https://doi.org/10.1002/mrm.28739.)
CO/computeVOP_classif_code
	VOP classification code used in computeVOP_CO and computeVOPi_CO
	nyc = 2; % not yet classified
	vop = 1; % classified as VOP
	dominated = 0; % classified as non-vop
CO/rQstar
	Solve the "generalized" eigenvalue problem max x^H.Q.x subject to max_Q*(x^H.Q*.x) <= 1
CO/spectralNorm
	Spectral norm computation across multiple SAR matrices using parfor
CO/CHtoRS
	Transforms a NxN complex Hermitian matrix into a 2Nx2N real symmetric matrix

CC/VOP_compression_iterative_VG
	VOP compression algorithm published in  (https://doi.org/10.1002/mrm.28739.) minimally modified for the puropose of the comparative study
	See https://sourceforge.net/projects/ relative-overestimation-vop/ for the orginal sources

CLU/VOPcompressionCLU
	Wrapper for VOPCompression.exe provided by Siemens Healthineers as part of Tim TX Step 2.3 product. 
	Please contact Siemens Healthineers to have access to VOPCompression.exe 

script_VOP
	Compute VOPs using the three compression methods computeVOPi_CO, VOP_compression_iterative_VG and VOPcompressionCLU and on the three test cases (Nova (8TX), Rapid (8TX) and Avanti2 (16 TX))
	To have access to the Avanti2 Q-matrices please mail to vincent.gras@cea.fr

script_VOPcheck
	run rQstar(VOP, Q) on each SAR matrix Q of the SAR model stored in the file passed as input

script_VOPcheckOverestimation
	Get SAR overestimation on each SAR matrix of the SAR model stored in the file passed as input
	
plotting_*
	Scripts to generate the figures
	!! some dependences are missing !! any request : mailto vincent.gras@cea.fr 

plot_ct
	plotting of computation time curves
	
plot_nVOP
	plotting number of VOPs curves

ct_*
	lolog analysis scripts
