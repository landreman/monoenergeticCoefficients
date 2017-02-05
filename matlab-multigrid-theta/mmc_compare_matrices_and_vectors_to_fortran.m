function mmc_compare_matrices_and_vectors_to_fortran(directory)

global N_levels levels fine_matrix rhs smoothing_option
global multigrid_restriction_matrices multigrid_prolongation_matrices

%addpath('~/petsc-3.6.3/share/petsc/matlab')
addpath('/Users/mattland/petsc/petsc-3.7.5/share/petsc/matlab')

% Compare restriction matrices:
for level = 1:(N_levels-1)
    filename = fullfile(directory,sprintf('mmc_restriction_matrix_level_%d.dat',level));
    if exist(filename,'file')==0
        fprintf('** Expected file %s does not exist.\n',filename)
    else
        matrix_fortran = PetscBinaryRead(filename);
        difference = full(max(max(abs(matrix_fortran - multigrid_restriction_matrices{level}))));
        if difference < 1e-8
            fprintf('  Restriction matrix on level %d agrees. (diff=%g)\n',level,difference)
        else
            fprintf('** Difference = %g between restriction matrices on level %d.',difference,level)
        end
    end
end


% Compare prolongation matrices:
for level = 1:(N_levels-1)
    filename = fullfile(directory,sprintf('mmc_prolongation_matrix_level_%d.dat',level));
    if exist(filename,'file')==0
        fprintf('** Expected file %s does not exist.\n',filename)
    else
        matrix_fortran = PetscBinaryRead(filename);
        difference = full(max(max(abs(matrix_fortran - multigrid_prolongation_matrices{level}))));
        if difference < 1e-8
            fprintf('  Prolongation matrix on level %d agrees. (diff=%g)\n',level,difference)
        else
            fprintf('** Difference = %g between prolongation matrices on level %d.',difference,level)
        end
    end
end

% Compare low-order matrices:
for level = 1:N_levels
    filename = fullfile(directory,sprintf('mmc_matrix_level_%d_whichMatrix_0.dat',level));
    if exist(filename,'file')==0
        fprintf('** Expected file %s does not exist.\n',filename)
    else
        matrix_fortran = PetscBinaryRead(filename);
        difference = full(max(max(abs(matrix_fortran - levels(level).low_order_matrix))));
        if difference < 1e-8
            fprintf('  Low order matrix on level %d agrees. (diff=%g)\n',level,difference)
        else
            fprintf('** Difference = %g between low order matrices on level %d.',difference,level)
        end
    end
end

% Compare high-order matrix on the fine level
filename = fullfile(directory,'mmc_matrix_level_1_whichMatrix_1.dat');
if exist(filename,'file')==0
    fprintf('** Expected file %s does not exist.\n',filename)
else
    matrix_fortran = PetscBinaryRead(filename);
    difference = full(max(max(abs(matrix_fortran - fine_matrix))));
    if difference < 1e-8
        fprintf('  High order matrix agrees. (diff=%g)\n',difference)
    else
        fprintf('** Difference = %g between high order matrices.',difference)
    end
end

% Compare right hand side vector
filename = fullfile(directory,'mmc_rhs.dat');
if exist(filename,'file')==0
    fprintf('** Expected file %s does not exist.\n',filename)
else
    matrix_fortran = PetscBinaryRead(filename);
    difference = full(max(max(abs(matrix_fortran - rhs))));
    if difference < 1e-8
        fprintf('  RHS vector agrees. (diff=%g)\n',difference)
    else
        fprintf('** Difference = %g between RHS vectors.',difference)
    end
end

% Compare smoothing matrices
for level = 1:(N_levels-1)
    switch smoothing_option
        case 1
            filename = fullfile(directory,sprintf('mmc_Jacobi_iteration_matrix_level_%d.dat',level));
            if exist(filename,'file')==0
                fprintf('** Expected file %s does not exist.\n',filename)
            else
                matrix_fortran = PetscBinaryRead(filename);
                difference = full(max(max(abs(matrix_fortran - levels(level).Jacobi_iteration_matrix))));
                if difference < 1e-8
                    fprintf('  Jacobi iteration matrix on level %d agrees. (diff=%g)\n',level,difference)
                else
                    fprintf('** Difference = %g between Jacobi iteration matrices on level %d.',difference,level)
                end
            end
        otherwise
            error('Invalid smoothing_option')
    end
end
return


% Check how many Jacobian matrices were saved:
for i=0:100
    filename = fullfile(directory,sprintf('sfincsBinary_iteration_%03d_whichMatrix_1',i));
    if exist(filename,'file')==0
        break
    end
    last_whichmatrix_1_filename = filename;
end
last_whichmatrix_1_found = i-1;
if i==0
    fprintf('No sfincsBinary_iteration_000_whichMatrix_1 file found.\n')
    return
end
fprintf('Iteration for last Jacobian matrix found: %d\n',last_whichmatrix_1_found)

% Check how many state vectors were saved:
% (So far this number has been identical to the number of matrices saved,
% but in principle it could be different if the SNES settings were
% different.)
for i=0:100
    filename = fullfile(directory,sprintf('sfincsBinary_iteration_%03d_stateVector',i));
    if exist(filename,'file')==0
        break
    end
    last_stateVector_filename = filename;
end
last_stateVector_found = i-1;
if i==0
    fprintf('No sfincsBinary_iteration_000_stateVector file found.\n')
    return
end
fprintf('Iteration for last state vector found: %d\n',last_stateVector_found)

% For the residual, always read iteration 000, since this initial residual
% is typically more interesting than later residuals.
filename = fullfile(directory,'sfincsBinary_iteration_000_residual');
fprintf('Attempting to read %s\n',filename)
residual_fortran = PetscBinaryRead(filename);

% For the preconditioner, always read iteration 000, since usually the
% preconditioner is not updated every iteration.
filename = fullfile(directory,'sfincsBinary_iteration_000_whichMatrix_0');
fprintf('Attempting to read %s\n',filename)
preconditionerMatrix_fortran = PetscBinaryRead(filename);

%filename = fullfile(directory,'sfincsBinary_iteration_000_whichMatrix_1');
filename = last_whichmatrix_1_filename;
fprintf('Attempting to read %s\n',filename)
Jacobian_fortran = PetscBinaryRead(filename);

%filename = fullfile(directory,'sfincsBinary_iteration_000_stateVector');
filename = last_stateVector_filename;
fprintf('Attempting to read %s\n',filename)
stateVector_fortran = PetscBinaryRead(filename);


figure(4)
clf
numRows = 2;
numCols = 2;

subplot(numRows,numCols,1)
plot(initialResidual,'.-','DisplayName','matlab')
hold on
plot(residual_fortran,'x:r','DisplayName','fortran')
legend show
title('Initial residual vector')

subplot(numRows,numCols,3)
plot(initialResidual - residual_fortran,'.-m')
title('Differences in initial residual vector')

subplot(numRows,numCols,2)
plot(stateVector,'.-','DisplayName','matlab')
hold on
plot(stateVector_fortran,'x:r','DisplayName','fortran')
legend show
title('Final stateVector')

subplot(numRows,numCols,4)
plot(stateVector - stateVector_fortran,'.-m')
title('Differences in final stateVector')

figure(5)
clf

numRows = 2;
numCols = 4;
th1 = 1e-8;
th2 = 1e-6;

subplot(numRows,numCols,1)
spy(abs(Jacobian)>th1)
title('Matlab Jacobian (last iteration)')
assignin('base','matlabJacobian_lastIteration',Jacobian)

subplot(numRows,numCols,2)
spy(abs(Jacobian_fortran)>th1)
title('Fortran Jacobian (last iteration)')
assignin('base','fortranJacobian_lastIteration',Jacobian_fortran)

differences = Jacobian_fortran-Jacobian;
assignin('base','JacobianDifferences',differences)

subplot(numRows,numCols,3)
%spy(abs(Jacobian_fortran-Jacobian)>th1)
spy(abs(differences)>1e-10,'k')
hold on
spy(abs(differences)>1e-9,'m')
spy(abs(differences)>1e-8,'b')
spy(abs(differences)>1e-7,'g')
spy(abs(differences)>1e-6,'r')
%title(['Differences > ',num2str(th2)])
title(['Differences > 1e-10,9,8,7,6'])
%title(['Differences > ',num2str(th1)])

subplot(numRows,numCols,4)
%spy(abs(Jacobian_fortran-Jacobian)>th2)
spy(abs(differences)>1e-5,'k')
hold on
spy(abs(differences)>1e-4,'m')
spy(abs(differences)>1e-3,'b')
spy(abs(differences)>1e-2,'g')
spy(abs(differences)>1e-1,'r')
%title(['Differences > ',num2str(th2)])
title(['Differences > 1e-5,4,3,2,1'])


%{
differencesMatter = abs(differences)>th2;
fprintf('Here come the differences:\n')
differences(differencesMatter)
%}

subplot(numRows,numCols,5)
spy(abs(preconditionerMatrix)>th1)
title('Matlab preconditioner (first iteration)')

subplot(numRows,numCols,6)
spy(abs(preconditionerMatrix_fortran)>th1)
title('Fortran preconditioner (first iteration)')

subplot(numRows,numCols,7)

%size(preconditionerMatrix_fortran)
%size(preconditionerMatrix)
spy(abs(preconditionerMatrix_fortran-preconditionerMatrix)>th1)
title(['Differences > ',num2str(th1)])

subplot(numRows,numCols,8)
spy(abs(preconditionerMatrix_fortran-preconditionerMatrix)>th2)
title(['Differences > ',num2str(th2)])

end