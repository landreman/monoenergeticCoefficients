global N_levels constraint_option restriction_option
global multigrid_restriction_matrices multigrid_prolongation_matrices
global levels Nperiods

multigrid_restriction_matrices = cell(N_levels-1,1);
multigrid_prolongation_matrices = cell(N_levels-1,1);

for level = 1:(N_levels-1)
    
    % ****************************************
    % Build prolongation operator
    % ****************************************
    xi_interpolation_matrix = sparse(mmc_nonperiodic_interpolation(levels(level+1).xi,levels(level).xi));
    theta_interpolation_matrix = sparse(mmc_periodic_interpolation(levels(level+1).theta,levels(level).theta,2*pi,2));
    zeta_interpolation_matrix = sparse(mmc_periodic_interpolation(levels(level+1).zeta,levels(level).zeta,2*pi/Nperiods,2));

    % Interpolation matrices should all have row sums of 1:
    assert(all(abs(sum(xi_interpolation_matrix,2)-1) < 1e-12))
    assert(all(abs(sum(theta_interpolation_matrix,2)-1) < 1e-12))
    assert(all(abs(sum(zeta_interpolation_matrix,2)-1) < 1e-12))
    
    %prolongation = kron(xi_interpolation_matrix,theta_interpolation_matrix);
    multigrid_prolongation_matrices{level} = kron(xi_interpolation_matrix,kron(theta_interpolation_matrix,zeta_interpolation_matrix));
    if constraint_option==1
        multigrid_prolongation_matrices{level}(end+1,end+1)=1;
    end
    
    % ****************************************
    % Build restriction operator
    % ****************************************
    restriction = multigrid_prolongation_matrices{level}';
    row_sums = sum(restriction,2);
    switch restriction_option
        case 1
            restriction = diag(1./row_sums) * restriction;
            % row sums should all be 1:
            assert(all(abs(sum(restriction,2)-1) < 1e-12))
        case 2
            restriction = restriction / row_sums(2); % Assume element 2 is in the interior.
        otherwise
            error('Invalid restriction_option')
    end
    multigrid_restriction_matrices{level} = restriction;
    
end
