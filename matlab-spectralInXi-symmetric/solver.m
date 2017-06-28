function [solution, totalNNZ] = solver(matrix, rhs, solutionMethod)

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
        
    otherwise
        error('Invalid value for solverType')
end

end