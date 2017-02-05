clear
%addpath('../src')

mmc_defaults()

% Change any parameters from the defaults here:

% For nu = 0.1,   you should get L11 = -0.0329,  L12 =  0.011
% For nu = 0.01,  you should get L11 = -0.0398,  L12 = -0.238
% For nu = 0.001, you should get L11 = -0.09,    L12 = -1.09

global nu E
nu = 0.01;
E = 0.03;

global Ntheta Nzeta Nxi
Ntheta = 15;
Nzeta = 15;
Nxi = 30;

global coarsen_theta coarsen_zeta coarsen_xi
coarsen_theta = false;
coarsen_zeta  = true;
coarsen_xi    = false;


mmc_main()


directory = ['/Users/mattland/monoenergeticCoefficients/fortran-multigrid-theta'];
%mmc_compareToFortran(fullfile(directory,'sfincsOutput.h5'))
mmc_compareMatricesAndVectorsToFortran(directory)
