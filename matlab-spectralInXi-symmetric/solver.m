function [solution, totalNNZ] = solver(matrix, preconditioner, rhs, solutionMethod, block_size)

% Options for solutionMethod:
% 1 = direct solver (backslash)
% 2 = direct solver (explicit LU)
% 3 = GMRES with no preconditioner

switch solutionMethod
    case 1
        % Direct solver
        fprintf('Using direct solver.\n')
        solution = matrix \ rhs;
        totalNNZ = -1; % Cannot be determined for this solutionMethod
        
    case 2
        fprintf('LU-factorizing the matrix.\n')
        [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(matrix);
        totalNNZ = nnz(preconditioner_L)+nnz(preconditioner_U);
        fprintf('nnz(L): %d,  nnz(U): %d,  total: %d\n',nnz(preconditioner_L), nnz(preconditioner_U), totalNNZ)
        solution = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * (rhs))));
        
        %[preconditioner_R, preconditioner_p, preconditioner_S] = chol(matrix);
        %totalNNZ = nnz(preconditioner_R);
        %fprintf('nnz(R): %d\n',totalNNZ)
        % Definition of permutation matrix S: if [R,p,S]=chol(A), then  R' R = S' A S
        %solution = preconditioner_S * (preconditioner_R \ ((preconditioner_R') \ (preconditioner_S' * rhs)));
        
    case 3
        % GMRES with no preconditioner
        totalNNZ = 0; % Memory required is negligible for this solutionMethod.
        restart = 100;
        tol = 1e-8;
        maxIterations = 200;
        
        fprintf('Beginning GMRES.\n')
        
        [solution,fl0,rr0,it0,rv0]=gmres(matrix,rhs,restart,tol,maxIterations/restart);
        switch fl0
            case 0
                fprintf('Converged!\n')
            case 1
                fprintf('Did not converge :(\n')
            case 2
                fprintf('Preconditioner was ill-conditioned\n')
            case 3
                fprintf('Stagnated :(\n')
        end
        figure(3)
        clf
        semilogy(rv0/rv0(1),'-o');
        xlabel('Iteration number');
        ylabel('Relative residual');
        title('Convergence of GMRES');
        
    case 4
        % GMRES with a preconditioner
        
        fprintf('LU-factorizing the preconditioner matrix.\n')
        [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditioner);
        totalNNZ = nnz(preconditioner_L)+nnz(preconditioner_U);
        fprintf('nnz(L): %d,  nnz(U): %d,  total: %d\n',nnz(preconditioner_L), nnz(preconditioner_U), totalNNZ)
        %solution = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * (rhs))));
        
        restart = 100;
        tol = 1e-8;
        maxIterations = 200;
        
        fprintf('Beginning GMRES.\n')
        x0 = zeros(size(rhs));
        [solution,fl0,rr0,it0,rv0]=gmres(matrix,rhs,restart,tol,maxIterations/restart,@apply_preconditioner, [], x0);
        switch fl0
            case 0
                fprintf('Converged!\n')
            case 1
                fprintf('Did not converge :(\n')
            case 2
                fprintf('Preconditioner was ill-conditioned\n')
            case 3
                fprintf('Stagnated :(\n')
        end
        figure(3)
        clf
        semilogy(rv0/rv0(1),'-o');
        xlabel('Iteration number');
        ylabel('Relative residual');
        title('Convergence of GMRES');
        
    case 5
        % MINRES with a preconditioner
        
        fprintf('LU-factorizing the preconditioner matrix.\n')
        [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditioner);
        totalNNZ = nnz(preconditioner_L)+nnz(preconditioner_U);
        fprintf('nnz(L): %d,  nnz(U): %d,  total: %d\n',nnz(preconditioner_L), nnz(preconditioner_U), totalNNZ)
        %solution = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * (rhs))));
        
        tol = 1e-8;
        maxIterations = 200;
        
        fprintf('Beginning MINRES.\n')
        [solution,fl0,rr0,it0,rv0] = minres(matrix,rhs,tol,maxIterations,@apply_preconditioner);
        switch fl0
            case 0
                fprintf('Converged!\n')
            case 1
                fprintf('Did not converge :(\n')
            case 2
                fprintf('Preconditioner was ill-conditioned\n')
            case 3
                fprintf('Stagnated :(\n')
            case 4
                fprintf('One of the scalar quantities calculated during minres became too small or too large to continue computing.\n')
            otherwise
                fprintf('Unexpected value for fl0: %g\n',fl0)
        end
        figure(3)
        clf
        semilogy(rv0/rv0(1),'-o');
        xlabel('Iteration number');
        ylabel('Relative residual');
        title('Convergence of MINRES');
        
    case 6
        % GMRES, preconditioning based on the 2x2 block structure
        
        matrix_a = preconditioner(1:block_size, 1:block_size);
        matrix_b = preconditioner(1:block_size, (block_size+1):end);
        matrix_c = preconditioner((block_size+1):end, 1:block_size);
        
        fprintf('LU-factorizing the preconditioner matrix.\n')
        [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditioner);
        totalNNZ = nnz(preconditioner_L)+nnz(preconditioner_U);
        fprintf('nnz(L): %d,  nnz(U): %d,  total: %d\n',nnz(preconditioner_L), nnz(preconditioner_U), totalNNZ)
        %solution = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * (rhs))));
        
        restart = 100;
        tol = 1e-8;
        maxIterations = 200;
        
        fprintf('Beginning GMRES.\n')
        x0 = zeros(size(rhs));
        [solution,fl0,rr0,it0,rv0]=gmres(matrix,rhs,restart,tol,maxIterations/restart,@apply_preconditioner, [], x0);
        switch fl0
            case 0
                fprintf('Converged!\n')
            case 1
                fprintf('Did not converge :(\n')
            case 2
                fprintf('Preconditioner was ill-conditioned\n')
            case 3
                fprintf('Stagnated :(\n')
        end
        figure(3)
        clf
        semilogy(rv0/rv0(1),'-o');
        xlabel('Iteration number');
        ylabel('Relative residual');
        title('Convergence of GMRES');
        
    otherwise
        error('Invalid value for solverType')
end

    function solnVector = apply_preconditioner(rhsVector)
        solnVector = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * rhsVector)));
    end
end