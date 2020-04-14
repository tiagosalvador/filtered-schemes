% This script compiles the mex files needed to run the solvers

addpath(genpath(pwd))
cd 'Solvers'/
mex -compatibleArrayDims loopMonotone.c
mex -compatibleArrayDims loopFilter.c
cd ..