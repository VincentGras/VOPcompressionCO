Iterative VOP algorithm.

Please cite 
Orzada S, Fielder TM, Quick HH, Ladd ME, 
"A local SAR compression algorithm with improved compression, speed and flexibility", 
https://Doi.org/10.1002/mrm.28739.

This new algorithm is an iterative expansion of the exisisting local SAR compression algorithm for MRI proposed by Lee et al. ( https://doi.org/10.1002/mrm.23140)

VOP_compression_Lee_SOR.m is an implementation of Lee's algorithm.

VOP_compression_iterative.m is the proposed new algorithm.

The number of VOPs is approximately halved with the new algorithm, while at the same time the compression time is reduced with speed-up factors of up to 2.5.

Code written by Stephan Orzada at German Cancer Research Center (DKFZ).