function solution = solver(matrix, rhs, solutionMethod)

% Options for solutionMethod:
% 1 = direct solver
% 2 = GMRES with no preconditioner

switch solutionMethod
    case 1
        % Direct solver
        fprintf('Using direct solver.\n')
        solution = matrix \ rhs;
        
    case 2
        % GMRES with no preconditioner
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