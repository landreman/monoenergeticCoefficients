function mmc_solver()

global levels N_levels solution N_smoothing Jacobi_omega
global gmres_restart gmres_maxit gmres_tol
global rhs fine_matrix smoothing_option
global coarse_L coarse_U coarse_P coarse_Q
global multigrid_restriction_matrices multigrid_prolongation_matrices

fprintf('Beginning GMRES.\n')

[solution,fl0,rr0,it0,rv0]=gmres(fine_matrix,rhs,gmres_restart,gmres_tol,gmres_maxit/gmres_restart,@apply_multigrid_cycle);
fprintf('After %d iterations, ',numel(rv0))
switch fl0
    case 0
        fprintf('converged!\n')
    case 1
        fprintf('did not converge :(\n')
    case 2
        fprintf('preconditioner was ill-conditioned\n')
    case 3
        fprintf('stagnated :(\n')
end

figure(7)
clf
semilogy(rv0/rv0(1),'-o');
xlabel('Iteration number');
ylabel('Relative residual');
title('Convergence of GMRES preconditioned with multigrid');

    function solution_vector = apply_multigrid_cycle(rhs_vector)
        fprintf(' Entering apply_multigrid_cycle.\n')

        levels(1).rhs_vector = rhs_vector;
        
        % Proceed from the fine level to the coarse level
        for level = 1:(N_levels-1)
            
            % Initialize smoothing for this level:
            switch smoothing_option
                case 1
                    levels(level).smoothing_shift = Jacobi_omega * (levels(level).Jacobi_D \ levels(level).rhs_vector);
                    levels(level).smoothing_shift(end) = 0;
                otherwise
                    error('Invalid smoothing_option')
            end
            levels(level).solution_vector = zeros(levels(level).matrixSize,1);
            
            % Pre-smoothing
            fprintf('  Pre-smoothing on level %d\n',level)
            switch smoothing_option
                case 1
                    for kk = 1:N_smoothing
                        levels(level).solution_vector = levels(level).Jacobi_iteration_matrix * levels(level).solution_vector + levels(level).smoothing_shift;
                    end
                otherwise
                    error('Invalid smoothing_option')
            end
            
            % Construct residual:
            levels(level).residual = levels(level).rhs_vector - levels(level).low_order_matrix * levels(level).solution_vector;
            
            % Restrict the residual to the next level:
            levels(level+1).rhs_vector = multigrid_restriction_matrices{level} * levels(level).residual;
        end
        
        fprintf('  Applying direct solver on the coarsest level.\n')
        levels(N_levels).solution_vector = coarse_Q * (coarse_U \ (coarse_L \ (coarse_P * levels(N_levels).rhs_vector)));
            
        % Proceed from the coarse level to the fine level
        for level = (N_levels-1):(-1):1

            % Update the fine-grid solution:
            levels(level).solution_vector = levels(level).solution_vector + multigrid_prolongation_matrices{level} * levels(level+1).solution_vector;
            
            % Pre-smoothing
            fprintf('  Post-smoothing on level %d\n',level)
            switch smoothing_option
                case 1
                    for kk = 1:N_smoothing
                        levels(level).solution_vector = levels(level).Jacobi_iteration_matrix * levels(level).solution_vector + levels(level).smoothing_shift;
                    end
                otherwise
                    error('Invalid smoothing_option')
            end
        end
        
        solution_vector = levels(1).solution_vector;
    end


end