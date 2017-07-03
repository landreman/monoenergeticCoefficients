filename = '/Users/mattland/monoenergeticCoefficients/fortran-spectralInXi-symmetric/mmc_rhs.dat';
rhs_fortran = PetscBinaryRead(filename);
rhs_matlab = problem.rhs;

filename = '/Users/mattland/monoenergeticCoefficients/fortran-spectralInXi-symmetric/mmc_solution.dat';
solution_fortran = PetscBinaryRead(filename);
solution_matlab = solution;

filename = '/Users/mattland/monoenergeticCoefficients/fortran-spectralInXi-symmetric/mmc_matrix.dat';
matrix_fortran = PetscBinaryRead(filename);
matrix_matlab = problem.matrix;

filename = '/Users/mattland/monoenergeticCoefficients/fortran-spectralInXi-symmetric/mmc_V.dat';
V_fortran = PetscBinaryRead(filename);
V_matlab = problem.V;

N=500;
V_fortran_small = full(V_fortran(1:N,1:N));
V_matlab_small  = full(V_matlab(1:N,1:N));

figure(1)
clf
plot(rhs_matlab,'.-','DisplayName','matlab')
hold on
plot(rhs_fortran,'.:r','DisplayName','fortran')
plot(rhs_matlab-rhs_fortran,'.:c','DisplayName','difference')
legend show
title('rhs')

fprintf('Max abs diff in rhs: %g\n',max(abs(rhs_matlab-rhs_fortran)))

figure(4)
clf
plot(solution_matlab,'.-','DisplayName','matlab')
hold on
plot(solution_fortran,'.:r','DisplayName','fortran')
plot(solution_matlab-solution_fortran,'.:c','DisplayName','difference')
legend show
title('solution')

fprintf('Max abs diff in solution: %g\n',max(abs(solution_matlab-solution_fortran)))

figure(2)
clf
numRows=1;
numCols=3;
th = 1e-10;

subplot(numRows,numCols,1)
spy(abs(matrix_matlab)>th)
title('matlab')

subplot(numRows,numCols,2)
spy(abs(matrix_fortran)>th)
title('fortran')

subplot(numRows,numCols,3)
spy(abs(matrix_matlab - matrix_fortran)>th)
title('difference')

figure(3)
clf
numRows=1;
numCols=3;
th = 1e-10;

subplot(numRows,numCols,1)
spy(abs(V_matlab)>th)
title('matlab')

subplot(numRows,numCols,2)
spy(abs(V_fortran)>th)
title('fortran')

subplot(numRows,numCols,3)
spy(abs(V_matlab - V_fortran)>th)
title('difference')

