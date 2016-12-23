clear

geometryParameters = struct(...
    'epsilon_t',-0.07053,...
    'epsilon_h',0.05067,...
    'iota',0.4542,...
    'G',3.7481,...
    'I',0,...
    'Nperiods',10,...
    'helicity_l',2);


%log10nus = linspace(-3,-1,3)
log10nus = (-5):0.2:(-1)
    
resolutionParameters = struct(...
    'Nalpha', 25,...
    'Nzeta', 29,...
    'Nxi', 70,...
    'includeConstraint', true);


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

log10nuConverged = log10nus(1);

quantitiesToRecord = {'flux','flow','time for solve','nnz(L)+nnz(U)','Num gmres iterations'};

linespecs = {'.-b','.-r','.-g','.-m','.-k','.-r','.:k','.:b','.-m'};

parametersToVary = {'log_{10}nu'};
abscissae = {log10nus};
convergeds = {log10nuConverged};

numQuantities = numel(quantitiesToRecord);
numParameters = numel(parametersToVary);
quantities = cell(numParameters,1);
quantities{1} = zeros(numel(log10nus), numQuantities);
parameterScanNum = 1;

numRuns = numel(log10nus);
runNum = 1;

% Vary nu, keeping other numerical parameters fixed.
for iii = 1:numel(log10nus)
    nu = 10 ^ log10nus(iii);
    fprintf('Beginning solve %d of %d. nu = %g\n',runNum,numRuns,nu); runNum = runNum+1;
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

figure(10)
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
