function mmc_create_grids()

global N_levels Nthetas Nzetas Nxis levels constraint_option Nperiods
global zetaMax xi_quadrature_option
global theta_derivative_option zeta_derivative_option xi_derivative_option
global preconditioner_theta_derivative_option preconditioner_zeta_derivative_option preconditioner_xi_derivative_option
global pitch_angle_scattering_option preconditioner_pitch_angle_scattering_option

levels = struct([]);
for level = 1:N_levels
    matrixSize = Nthetas(level)*Nzetas(level)*Nxis(level);
    if constraint_option == 1
        matrixSize = matrixSize + 1;
    end
    levels(level).matrixSize = matrixSize;
    
    % *************************************************************
    % Build theta grid for main matrix
    % *************************************************************
    
    switch theta_derivative_option
        case 2
            fprintf('df/dtheta derivative: centered differences, 1 point on each side.\n')
            derivative_option_plus = 0;
            derivative_option_minus = derivative_option_plus;
        case 3
            fprintf('df/dtheta derivative: centered differences, 2 points on each side.\n')
            derivative_option_plus = 10;
            derivative_option_minus = derivative_option_plus;
        case 4
            fprintf('df/dtheta derivative: upwinded differences, 0 points on one side, 1 point on the other side.\n')
            derivative_option_plus  = 30;
            derivative_option_minus = 40;
        case 5
            fprintf('df/dtheta derivative: upwinded differences, 0 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 50;
            derivative_option_minus = 60;
        case 6
            fprintf('df/dtheta derivative: upwinded differences, 1 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 80;
            derivative_option_minus = 90;
        case 7
            fprintf('df/dtheta derivative: upwinded differences, 1 point on one side, 3 points on the other side.\n')
            derivative_option_plus  = 100;
            derivative_option_minus = 110;
        case 8
            fprintf('df/dtheta derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
            derivative_option_plus  = 120;
            derivative_option_minus = 130;
        otherwise
            error('Invalid theta_derivative_option')
    end
    quadrature_option = 0;
    [levels(level).theta, levels(level).thetaWeights, levels(level).ddtheta_plus, ~]  = mmc_uniformDiffMatrices(Nthetas(level), ...
        0, 2*pi, derivative_option_plus, quadrature_option);
    [~, ~, levels(level).ddtheta_minus, ~] = mmc_uniformDiffMatrices(Nthetas(level), ...
        0, 2*pi, derivative_option_minus, quadrature_option);
    
    % *************************************************************
    % Build theta grid for preconditioner
    % *************************************************************
    
    call_uniform_diff_matrices = true;
    switch abs(preconditioner_theta_derivative_option)
        case 0
            fprintf('Preconditioner df/dtheta derivative is dropped.\n')
            call_uniform_diff_matrices = false;
            ddtheta_plus = zeros(Nthetas(level));
            ddtheta_minus = zeros(Nthetas(level));
        case 2
            fprintf('Preconditioner df/dtheta derivative: centered differences, 1 point on each side.\n')
            derivative_option_plus = 0;
            derivative_option_minus = derivative_option_plus;
        case 3
            fprintf('Preconditioner df/dtheta derivative: centered differences, 2 points on each side.\n')
            derivative_option_plus = 10;
            derivative_option_minus = derivative_option_plus;
        case 4
            fprintf('Preconditioner df/dtheta derivative: upwinded differences, 0 points on one side, 1 point on the other side.\n')
            derivative_option_plus  = 30;
            derivative_option_minus = 40;
        case 5
            fprintf('Preconditioner df/dtheta derivative: upwinded differences, 0 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 50;
            derivative_option_minus = 60;
        case 6
            fprintf('Preconditioner df/dtheta derivative: upwinded differences, 1 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 80;
            derivative_option_minus = 90;
        case 7
            fprintf('Preconditioner df/dtheta derivative: upwinded differences, 1 point on one side, 3 points on the other side.\n')
            derivative_option_plus  = 100;
            derivative_option_minus = 110;
        case 8
            fprintf('Preconditioner df/dtheta derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
            derivative_option_plus  = 120;
            derivative_option_minus = 130;
        otherwise
            error('Invalid preconditioner_theta_derivative_option')
    end
    if call_uniform_diff_matrices
        quadrature_option = 0;
        [~, ~, ddtheta_plus, ~]  = mmc_uniformDiffMatrices(Nthetas(level), ...
            0, 2*pi, derivative_option_plus, quadrature_option);
        [~, ~, ddtheta_minus, ~] = mmc_uniformDiffMatrices(Nthetas(level), ...
            0, 2*pi, derivative_option_minus, quadrature_option);
    end
    if preconditioner_theta_derivative_option<0
        ddtheta_plus = diag(diag(ddtheta_plus));
        ddtheta_minus = diag(diag(ddtheta_minus));
        fprintf('  But only the diagonal is kept.\n')
    end
    levels(level).ddtheta_plus_preconditioner = ddtheta_plus;
    levels(level).ddtheta_minus_preconditioner = ddtheta_minus;
    
    
    % *************************************************************
    % Build zeta grid for main matrix
    % *************************************************************
    
    zetaMax = 2*pi/Nperiods;
    
    switch zeta_derivative_option
        case 2
            fprintf('df/dzeta derivative: centered differences, 1 point on each side.\n')
            derivative_option_plus = 0;
            derivative_option_minus = derivative_option_plus;
        case 3
            fprintf('df/dzeta derivative: centered differences, 2 points on each side.\n')
            derivative_option_plus = 10;
            derivative_option_minus = derivative_option_plus;
        case 4
            fprintf('df/dzeta derivative: upwinded differences, 0 points on one side, 1 point on the other side.\n')
            derivative_option_plus  = 30;
            derivative_option_minus = 40;
        case 5
            fprintf('df/dzeta derivative: upwinded differences, 0 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 50;
            derivative_option_minus = 60;
        case 6
            fprintf('df/dzeta derivative: upwinded differences, 1 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 80;
            derivative_option_minus = 90;
        case 7
            fprintf('df/dzeta derivative: upwinded differences, 1 point on one side, 3 points on the other side.\n')
            derivative_option_plus  = 100;
            derivative_option_minus = 110;
        case 8
            fprintf('df/dzeta derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
            derivative_option_plus  = 120;
            derivative_option_minus = 130;
        otherwise
            error('Invalid zeta_derivative_option')
    end
    quadrature_option = 0;
    [levels(level).zeta, levels(level).zetaWeights, levels(level).ddzeta_plus, ~]  = mmc_uniformDiffMatrices(Nzetas(level), 0, zetaMax, derivative_option_plus, quadrature_option);
    [~   , ~, levels(level).ddzeta_minus, ~] = mmc_uniformDiffMatrices(Nzetas(level), 0, zetaMax, derivative_option_minus, quadrature_option);
    
    levels(level).zetaWeights = levels(level).zetaWeights * Nperiods;
    assert(abs(sum(levels(level).zetaWeights)-2*pi) < 1e-12)
    
    
    % *************************************************************
    % Build zeta grid for preconditioner
    % *************************************************************
    
    switch abs(preconditioner_zeta_derivative_option)
        case 2
            fprintf('Preconditioner df/dzeta derivative: centered differences, 1 point on each side.\n')
            derivative_option_plus = 0;
            derivative_option_minus = derivative_option_plus;
        case 3
            fprintf('Preconditioner df/dzeta derivative: centered differences, 2 points on each side.\n')
            derivative_option_plus = 10;
            derivative_option_minus = derivative_option_plus;
        case 4
            fprintf('Preconditioner df/dzeta derivative: upwinded differences, 0 points on one side, 1 point on the other side.\n')
            derivative_option_plus  = 30;
            derivative_option_minus = 40;
        case 5
            fprintf('Preconditioner df/dzeta derivative: upwinded differences, 0 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 50;
            derivative_option_minus = 60;
        case 6
            fprintf('Preconditioner df/dzeta derivative: upwinded differences, 1 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 80;
            derivative_option_minus = 90;
        case 7
            fprintf('Preconditioner df/dzeta derivative: upwinded differences, 1 point on one side, 3 points on the other side.\n')
            derivative_option_plus  = 100;
            derivative_option_minus = 110;
        case 8
            fprintf('Preconditioner df/dzeta derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
            derivative_option_plus  = 120;
            derivative_option_minus = 130;
        otherwise
            error('Invalid preconditioner_zeta_derivative_option')
    end
    quadrature_option = 0;
    [zeta, zetaWeights, ddzeta_plus, ~]  = mmc_uniformDiffMatrices(Nzetas(level), 0, zetaMax, derivative_option_plus, quadrature_option);
    [~   , ~, ddzeta_minus, ~] = mmc_uniformDiffMatrices(Nzetas(level), 0, zetaMax, derivative_option_minus, quadrature_option);
    
    if preconditioner_zeta_derivative_option<0
        fprintf('  But only the diagonal is kept.\n')
        ddzeta_plus = diag(diag(ddzeta_plus));
        ddzeta_minus = diag(diag(ddzeta_minus));
    end
    levels(level).ddzeta_plus_preconditioner = ddzeta_plus;
    levels(level).ddzeta_minus_preconditioner = ddzeta_minus;
    
    % *************************************************************
    % Initialize the arrays for the magnetic field strength B:
    % *************************************************************
    
    [levels(level).zeta2D, levels(level).theta2D] = meshgrid(levels(level).zeta,levels(level).theta);
    [levels(level).B, levels(level).dBdtheta, levels(level).dBdzeta] = mmc_geometry(levels(level).theta2D, levels(level).zeta2D);
    
    
    
    % *************************************************************
    % Generate xi integration weights. We'll deal with the differentiation
    % matrices later.
    % *************************************************************
    
    derivative_option=12;
    [levels(level).xi, levels(level).xiWeights, ~, ~]  = mmc_uniformDiffMatrices(Nxis(level), -1, 1, derivative_option, xi_quadrature_option);
    
    % *************************************************************
    % Build xi grid for main matrix
    % *************************************************************
    
    switch xi_derivative_option
        case 2
            fprintf('df/dxi derivative: centered differences, 1 point on each side.\n')
            derivative_option_plus = 2;
            derivative_option_minus = derivative_option_plus;
        case 3
            fprintf('df/dxi derivative: centered differences, 2 points on each side.\n')
            derivative_option_plus = 12;
            derivative_option_minus = derivative_option_plus;
        case 4
            fprintf('df/dxi derivative: upwinded differences, 0 points on one side, 1 point on the other side.\n')
            derivative_option_plus  = 32;
            derivative_option_minus = 42;
        case 5
            fprintf('df/dxi derivative: upwinded differences, 0 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 52;
            derivative_option_minus = 62;
        case 6
            fprintf('df/dxi derivative: upwinded differences, 1 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 82;
            derivative_option_minus = 92;
        case 7
            fprintf('df/dxi derivative: upwinded differences, 1 point on one side, 3 points on the other side.\n')
            derivative_option_plus  = 102;
            derivative_option_minus = 112;
        case 8
            fprintf('df/dxi derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
            fprintf('  High accuracy at the ends but imperfect upwinding.\n')
            derivative_option_plus  = 122;
            derivative_option_minus = 132;
        case 9
            fprintf('df/dxi derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
            fprintf('  Upwinding all the way to the ends but nonideal accuracy.\n')
            derivative_option_plus  = 123;
            derivative_option_minus = 133;
        otherwise
            error('Invalid xi_derivative_option')
    end
    quadrature_option = xi_quadrature_option;
    [~, ~, levels(level).ddxi_plus, ~]  = mmc_uniformDiffMatrices(Nxis(level), -1, 1, derivative_option_plus,  quadrature_option);
    [~, ~, levels(level).ddxi_minus, ~] = mmc_uniformDiffMatrices(Nxis(level), -1, 1, derivative_option_minus, quadrature_option);
    
    % *************************************************************
    % Build xi grid for preconditioner
    % *************************************************************
    
    call_uniform_diff_matrices = true;
    switch abs(preconditioner_xi_derivative_option)
        case 0
            fprintf('Preconditioner df/dxi derivative is dropped.\n')
            call_uniform_diff_matrices = false;
            ddxi_plus  = zeros(Nxis(level));
            ddxi_minus = zeros(Nxis(level));
        case 2
            fprintf('Preconditioner df/dxi derivative: centered differences, 1 point on each side.\n')
            derivative_option_plus = 2;
            derivative_option_minus = derivative_option_plus;
        case 3
            fprintf('Preconditioner df/dxi derivative: centered differences, 2 points on each side.\n')
            derivative_option_plus = 12;
            derivative_option_minus = derivative_option_plus;
        case 4
            fprintf('Preconditioner df/dxi derivative: upwinded differences, 0 points on one side, 1 point on the other side.\n')
            derivative_option_plus  = 32;
            derivative_option_minus = 42;
        case 5
            fprintf('Preconditioner df/dxi derivative: upwinded differences, 0 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 52;
            derivative_option_minus = 62;
        case 6
            fprintf('Preconditioner df/dxi derivative: upwinded differences, 1 points on one side, 2 points on the other side.\n')
            derivative_option_plus  = 82;
            derivative_option_minus = 92;
        case 7
            fprintf('Preconditioner df/dxi derivative: upwinded differences, 1 point on one side, 3 points on the other side.\n')
            derivative_option_plus  = 102;
            derivative_option_minus = 112;
        case 8
            fprintf('Preconditioner df/dxi derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
            fprintf('  High accuracy at the ends but imperfect upwinding.\n')
            derivative_option_plus  = 122;
            derivative_option_minus = 132;
        case 9
            fprintf('Preconditioner df/dxi derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
            fprintf('  Upwinding all the way to the ends but nonideal accuracy.\n')
            derivative_option_plus  = 123;
            derivative_option_minus = 133;
        otherwise
            error('Invalid preconditioner_xi_derivative_option')
    end
    if call_uniform_diff_matrices
        quadrature_option = xi_quadrature_option;
        [~, ~, ddxi_plus, ~]  = mmc_uniformDiffMatrices(Nxis(level), -1, 1, derivative_option_plus,  quadrature_option);
        [~, ~, ddxi_minus, ~] = mmc_uniformDiffMatrices(Nxis(level), -1, 1, derivative_option_minus, quadrature_option);
    end
    if preconditioner_xi_derivative_option<0
        ddxi_plus = diag(diag(ddxi_plus));
        ddxi_minus = diag(diag(ddxi_minus));
        fprintf('  But only the diagonal is kept.\n')
    end
    levels(level).ddxi_plus_preconditioner = ddxi_plus;
    levels(level).ddxi_minus_preconditioner = ddxi_minus;
    
    
    
    % *************************************************************
    % Build pitch angle scattering matrix for main matrix
    % *************************************************************

    switch pitch_angle_scattering_option
        case 2
            fprintf('Pitch angle scattering operator: centered differences, 1 point on each side.\n')
            derivative_option = 2;
        case 3
            fprintf('Pitch angle scattering operator: centered differences, 2 points on each side.\n')
            derivative_option = 12;
        otherwise
            error('Invalid pitch_angle_scattering_option')
    end
    [~, ~, ddxi, d2dxi2]  = mmc_uniformDiffMatrices(Nxis(level), -1, 1, derivative_option,  quadrature_option);
    xi = levels(level).xi;
    levels(level).pitch_angle_scattering_operator = 0.5*diag(1-xi.^2)*d2dxi2 - diag(xi)*ddxi;
    
    
    % *************************************************************
    % Build pitch angle scattering matrix for preconditioner
    % *************************************************************

    call_uniform_diff_matrices = true;
    switch abs(preconditioner_pitch_angle_scattering_option)
        case 0
            fprintf('Pitch angle scattering operator is dropped in the preconditioner.\n')
            call_uniform_diff_matrices = false;
            levels(level).pitch_angle_scattering_operator_preconditioner  = zeros(Nxis(level));
        case 2
            fprintf('Preconditioner pitch angle scattering operator: centered differences, 1 point on each side.\n')
            derivative_option = 2;
        case 3
            fprintf('Preconditioner pitch angle scattering operator: centered differences, 2 points on each side.\n')
            derivative_option = 12;
        otherwise
            error('Invalid preconditioner_pitch_angle_scattering_option')
    end
    if call_uniform_diff_matrices
        [~, ~, ddxi, d2dxi2]  = mmc_uniformDiffMatrices(Nxis(level), -1, 1, derivative_option,  quadrature_option);
        xi = levels(level).xi;
        pitch_angle_scattering_operator = 0.5*diag(1-xi.^2)*d2dxi2 - diag(xi)*ddxi;
        if preconditioner_pitch_angle_scattering_option<0
            pitch_angle_scattering_operator = diag(diag(pitch_angle_scattering_operator));
            fprintf('  But only the diagonal is kept.\n')
        end
        levels(level).pitch_angle_scattering_operator_preconditioner = pitch_angle_scattering_operator;
    end
    
end

global VPrime FSAB2
VPrime = levels(1).thetaWeights' * (1./levels(1).B.^2) * levels(1).zetaWeights;
FSAB2 = (levels(1).thetaWeights' * ones(size(levels(1).B)) * levels(1).zetaWeights) / VPrime;

end
