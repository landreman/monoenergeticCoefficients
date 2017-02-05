global N_levels levels fine_matrix smoothing_option Jacobi_omega
global flux flow

mmc_set_grid_resolutions()

mmc_create_grids()

mmc_restriction_prolongation_matrices()

% Build the high-order matrix on the fine level
which_matrix = 1;
level = 1;
fine_matrix = mmc_populate_matrix(which_matrix, level);

% Build the low-order matrices on every level
which_matrix = 0;
for level = 1:N_levels
    levels(level).low_order_matrix = mmc_populate_matrix(which_matrix, level);
end

% Build matrices used for smoothing:
for level = 1:(N_levels-1)
    switch smoothing_option
        case 1
            % Jacobi
            
            matrixSize = levels(level).matrixSize;
            D = diag(diag(levels(level).low_order_matrix));
            D(end,end)=1;
            levels(level).Jacobi_D = D;
            U_plus_L = levels(level).low_order_matrix - D;
            %Jacobi_iteration_matrix = (1-omega)*speye(matrixSize) - omega*(D \ U_plus_L);
            Jacobi_iteration_matrix = sparse((1-Jacobi_omega)*speye(matrixSize) - Jacobi_omega*(D \ U_plus_L));
            
            Jacobi_iteration_matrix(end,:)=0;
            Jacobi_iteration_matrix(:,end)=0;
            
            levels(level).Jacobi_iteration_matrix = Jacobi_iteration_matrix;
            
        case 2
            % Line smoothing in zeta
        otherwise
            error('Invalid smoothing_option')
    end
end

% LU factorize the low-order matrix on the coarsest level:
tic
fprintf('LU-factorizing the low-order matrix on the coarsest level.\n')
global coarse_L coarse_U coarse_P coarse_Q
[coarse_L, coarse_U, coarse_P, coarse_Q] = lu(levels(N_levels).low_order_matrix);
totalNNZ = nnz(coarse_L)+nnz(coarse_U);
fprintf('nnz(L): %d,  nnz(U): %d,  total: %d\n',nnz(coarse_L), nnz(coarse_U), totalNNZ)
fprintf('Took %g sec.\n',toc)

mmc_populate_RHS()

mmc_solver()

mmc_diagnostics()

fprintf('Diagnostics: flux (L11) = %g, flow (L21) = %g.\n',flux,flow);

