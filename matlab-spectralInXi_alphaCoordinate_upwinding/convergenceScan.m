clear

geometryParameters = struct(...
    'epsilon_t',-0.07053,...
    'epsilon_h',0.05067,...
    'iota',0.4542,...
    'G',3.7481,...
    'I',0,...
    'Nperiods',10,...
    'helicity_l',2);

%{
nu=0.1;

NalphaConverged = 15;
Nalphas = 11:2:33;

NzetaConverged = 15;
Nzetas = 11:2:33;

NxiConverged = 15;
Nxis = 12:30;
%Nxis = round(linspace(30,100,15));
%Nxis = round(linspace(10,30,15));
%}

%{
nu=0.01;
    
NalphaConverged = 15;
Nalphas = 11:2:35;

NzetaConverged = 15;
Nzetas = 11:2:35;

NxiConverged = 30;
Nxis = 25:60;
%Nxis = round(linspace(30,100,15));
%Nxis = round(linspace(10,30,15));
%}

nu=0.001;
    
NalphaConverged = 25;
Nalphas = 25:2:51;

NzetaConverged = 29;
Nzetas = 21:2:59;

NxiConverged = 50;
%Nxis = 25:60;
Nxis = round(linspace(35,100,15));
%Nxis = round(linspace(10,30,15));


%{
% Collisionality, using the same normalization as in SFINCS:
nuPrime = 0.1;

% For comparison with Hakan's convention, scale the collisionality by nu_D(1):
Psi_Chandra = (erf(1) - 2/sqrt(pi)*exp(-1)) / 2;
nuD = 3*sqrt(pi)/4*(erf(1) - Psi_Chandra);
nu = nuPrime * nuD;
%}
%nu=0.01;


discretizationParameters = struct(...
    'zeta_derivative_option', 8, ...
    'alpha_interpolation_stencil', 4, ...
    'include_xi_pentadiagonal_terms', true);

preconditionerDiscretizationParameters = struct(...
    'zeta_derivative_option', 4, ...
    'alpha_interpolation_stencil', 2, ...
    'include_xi_pentadiagonal_terms', true);

includeConstraint = true;

solutionMethod = 4;
% 1 = sparse direct solver (backslash)
% 2 = sparse direct solver (explicit LU decomposition, so memory can be monitored.)
% 3 = GMRES with no preconditioning
% 4 = GMRES with preconditioning

% **********************************************************
% End of input parameters.
% **********************************************************

switch discretizationParameters.zeta_derivative_option
    case {2,4}
        discretizationParameters.buffer_zeta_points_on_each_side = 1;
    case {3,5,6}
        discretizationParameters.buffer_zeta_points_on_each_side = 2;
    case {7,8}
        discretizationParameters.buffer_zeta_points_on_each_side = 3;
    otherwise
        error('Invalid discretizationParameters.zeta_derivative_option')
end
preconditionerDiscretizationParameters.buffer_zeta_points_on_each_side = discretizationParameters.buffer_zeta_points_on_each_side;



quantitiesToRecord = {'flux','flow','time for LU','nnz(L)+nnz(U)','Num gmres iterations'};

linespecs = {'.-b','.-r','.-g','.-m','.-k','.-r','.:k','.:b','.-m'};

parametersToVary = {'N\alpha','N\zeta','N\xi'};
abscissae = {Nalphas, Nzetas, Nxis};
convergeds = {NalphaConverged, NzetaConverged, NxiConverged};

numQuantities = numel(quantitiesToRecord);
numParameters = numel(parametersToVary);
quantities = cell(numParameters,1);
quantities{1} = zeros(numel(Nalphas), numQuantities);
quantities{2} = zeros(numel(Nzetas), numQuantities);
quantities{3} = zeros(numel(Nxis), numQuantities);
parameterScanNum = 1;

numRuns = numel(Nalphas)+numel(Nzetas)+numel(Nxis);
runNum = 1;

% Vary Nalpha, keeping other numerical parameters fixed.
for iii = 1:numel(Nalphas)
    resolutionParameters = struct(...
        'Nalpha', Nalphas(iii),...
        'Nzeta', NzetaConverged,...
        'Nxi', NxiConverged,...
        'includeConstraint', includeConstraint);
    fprintf('Beginning solve %d of %d.\n',runNum,numRuns); runNum = runNum+1;
    problem = assembleMatrix(resolutionParameters, nu, geometryParameters, discretizationParameters);
    preconditioner = assembleMatrix(resolutionParameters, nu, geometryParameters, preconditionerDiscretizationParameters);
    tic
    [solution, totalNNZ, num_iterations] = solver(problem.matrix, problem.rhs, preconditioner.matrix, solutionMethod);
    outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
    quantities{parameterScanNum}(iii,1)=outputs.flux;
    quantities{parameterScanNum}(iii,2)=outputs.flow;
    quantities{parameterScanNum}(iii,3)=toc;
    quantities{parameterScanNum}(iii,4)=totalNNZ;
    quantities{parameterScanNum}(iii,5)=num_iterations;
end
parameterScanNum = parameterScanNum+1;

% Vary Nzeta, keeping other numerical parameters fixed.
for iii = 1:numel(Nzetas)
    resolutionParameters = struct(...
        'Nalpha', NalphaConverged,...
        'Nzeta', Nzetas(iii),...
        'Nxi', NxiConverged,...
        'includeConstraint', includeConstraint);
    fprintf('Beginning solve %d of %d.\n',runNum,numRuns); runNum = runNum+1;
    problem = assembleMatrix(resolutionParameters, nu, geometryParameters, discretizationParameters);
    preconditioner = assembleMatrix(resolutionParameters, nu, geometryParameters, preconditionerDiscretizationParameters);
    tic
    [solution, totalNNZ, num_iterations] = solver(problem.matrix, problem.rhs, preconditioner.matrix, solutionMethod);
    outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
    quantities{parameterScanNum}(iii,1)=outputs.flux;
    quantities{parameterScanNum}(iii,2)=outputs.flow;
    quantities{parameterScanNum}(iii,3)=toc;
    quantities{parameterScanNum}(iii,4)=totalNNZ;
    quantities{parameterScanNum}(iii,5)=num_iterations;
end
parameterScanNum = parameterScanNum+1;

% Vary Nxi, keeping other numerical parameters fixed.
for iii = 1:numel(Nxis)
    resolutionParameters = struct(...
        'Nalpha', NalphaConverged,...
        'Nzeta', NzetaConverged,...
        'Nxi', Nxis(iii),...
        'includeConstraint', includeConstraint);
    fprintf('Beginning solve %d of %d.\n',runNum,numRuns); runNum = runNum+1;
    problem = assembleMatrix(resolutionParameters, nu, geometryParameters, discretizationParameters);
    preconditioner = assembleMatrix(resolutionParameters, nu, geometryParameters, preconditionerDiscretizationParameters);
    tic
    [solution, totalNNZ, num_iterations] = solver(problem.matrix, problem.rhs, preconditioner.matrix, solutionMethod);
    outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
    quantities{parameterScanNum}(iii,1)=outputs.flux;
    quantities{parameterScanNum}(iii,2)=outputs.flow;
    quantities{parameterScanNum}(iii,3)=toc;
    quantities{parameterScanNum}(iii,4)=totalNNZ;
    quantities{parameterScanNum}(iii,5)=num_iterations;
end
parameterScanNum = parameterScanNum+1;

maxs=ones(numQuantities,1)*(-Inf);
mins=ones(numQuantities,1)*(Inf);
for iParameter = 1:numParameters
    maxs = max([maxs, quantities{iParameter}(:,:)'],[],2);
    mins = min([mins, quantities{iParameter}(:,:)'],[],2);
end

figure(9)
clf
numRows = numQuantities;
numCols = numParameters;
for iQuantity = 1:numQuantities
    if maxs(iQuantity) <= mins(iQuantity)
        maxs(iQuantity) = mins(iQuantity)+1;
    end

    for iParameter = 1:numParameters
        subplot(numRows, numCols, iParameter  + (iQuantity - 1)*numParameters)
        plot(abscissae{iParameter}, quantities{iParameter}(:,iQuantity)', linespecs{iQuantity})
        hold on
        plot([convergeds{iParameter}, convergeds{iParameter}], [mins(iQuantity),maxs(iQuantity)],'k')
        ylim([mins(iQuantity), maxs(iQuantity)])
        xlabel(parametersToVary{iParameter})
        ylabel(quantitiesToRecord{iQuantity})
    end
end
