function [solution, totalNNZ, num_iterations] = solver(matrix, rhs, preconditioner, solutionMethod)

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
        num_iterations = 1;
        
    case 2
        fprintf('LU-factorizing the matrix.\n')
        [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(matrix);
        totalNNZ = nnz(preconditioner_L)+nnz(preconditioner_U);
        fprintf('nnz(L): %d,  nnz(U): %d,  total: %d\n',nnz(preconditioner_L), nnz(preconditioner_U), totalNNZ)
        solution = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * (rhs))));
        num_iterations = 1;
        
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
        
        num_iterations = numel(rv0);
    case 4
        % GMRES with preconditioner
        
        fprintf('LU-factorizing the preconditioner.\n')
        [preconditioner_L, preconditioner_U, preconditioner_P, preconditioner_Q] = lu(preconditioner);
        totalNNZ = nnz(preconditioner_L)+nnz(preconditioner_U);
        fprintf('nnz(L): %d,  nnz(U): %d,  total: %d\n',nnz(preconditioner_L), nnz(preconditioner_U), totalNNZ)

        restart = 3;
        tol = 1e-8;
        maxIterations = 200;
        
        fprintf('Beginning GMRES.\n')
        
        [solution,fl0,rr0,it0,rv0]=gmres(matrix,rhs,restart,tol,maxIterations/restart,@applyPreconditioner);
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
        
        num_iterations = numel(rv0);

    otherwise
        error('Invalid value for solverType')
end

    function solnVector = applyPreconditioner(rhsVector)
        solnVector = preconditioner_Q * (preconditioner_U \ (preconditioner_L \ (preconditioner_P * rhsVector)));
    end

end