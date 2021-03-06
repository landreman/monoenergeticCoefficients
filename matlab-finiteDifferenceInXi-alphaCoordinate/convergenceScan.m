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
NalphaConverged = 15;
Nalphas = 11:2:39;

NzetaConverged = 19;
Nzetas = 15:2:33;

NxiConverged = 50;
Nxis = 25:60;
%Nxis = round(linspace(30,100,15));
%Nxis = round(linspace(10,30,15));
%}
%{
NalphaConverged = 15;
Nalphas = 11:2:39;

NzetaConverged = 15;
Nzetas = 11:2:39;

NxiConverged = 50;
Nxis = 25:60;
%Nxis = round(linspace(30,100,15));
%Nxis = round(linspace(10,30,15));
%}
%{
NalphaConverged = 15;
Nalphas = 11:2:31;

NzetaConverged = 15;
Nzetas = 11:2:31;

NxiConverged = 16;
Nxis = 15:30;
%}
%{
NalphaConverged = 15;
Nalphas = 11:2:31;

NzetaConverged = 15;
Nzetas = 11:2:31;

NxiConverged = 30;
Nxis = 25:60;
%}

NalphaConverged = 19;
Nalphas = 11:2:39;

NzetaConverged = 19;
Nzetas = 11:2:39;

NxiConverged = 45;
Nxis = round(linspace(30,90,20));

    
    %{
% Collisionality, using the same normalization as in SFINCS:
nuPrime = 0.1;

% For comparison with Hakan's convention, scale the collisionality by nu_D(1):
Psi_Chandra = (erf(1) - 2/sqrt(pi)*exp(-1)) / 2;
nuD = 3*sqrt(pi)/4*(erf(1) - Psi_Chandra);
nu = nuPrime * nuD;
%}
nu=0.001;
E = 0;

ExB_alpha_option_1 = 0;
zeta_preconditioner_option_1  = 6;
xi_preconditioner_option_1    = 0;

ExB_alpha_option_2 = 0;
zeta_preconditioner_option_2  = 2;
xi_preconditioner_option_2    = 2;

includeConstraint = true;
solutionMethod = -4;

quantitiesToRecord = {'flux','flow','time for LU','nnz(L)+nnz(U)'};

linespecs = {'.-b','.-r','.-g','.-m','.:c','.-r','.:k','.:b','.-m'};

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
    problem = assembleMatrix(resolutionParameters, nu, E, geometryParameters, ...
        ExB_alpha_option_1,zeta_preconditioner_option_1,xi_preconditioner_option_1);
    if abs(solutionMethod)>2
        preconditioner = assembleMatrix(resolutionParameters, nu, E, geometryParameters, ...
            ExB_alpha_option_2,zeta_preconditioner_option_2,xi_preconditioner_option_2);
    else
        preconditioner=problem;
    end
    tic
    [solution, totalNNZ] = solver(problem.matrix, problem.rhs, preconditioner.matrix, solutionMethod);
    outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
    quantities{parameterScanNum}(iii,1)=outputs.flux;
    quantities{parameterScanNum}(iii,2)=outputs.flow;
    quantities{parameterScanNum}(iii,3)=toc;
    quantities{parameterScanNum}(iii,4)=totalNNZ;
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
    problem = assembleMatrix(resolutionParameters, nu, E, geometryParameters, ...
        ExB_alpha_option_1,zeta_preconditioner_option_1,xi_preconditioner_option_1);
    if abs(solutionMethod)>2
        preconditioner = assembleMatrix(resolutionParameters, nu, E, geometryParameters, ...
            ExB_alpha_option_2,zeta_preconditioner_option_2,xi_preconditioner_option_2);
    else
        preconditioner=problem;
    end
    tic
    [solution, totalNNZ] = solver(problem.matrix, problem.rhs, preconditioner.matrix, solutionMethod);
    outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
    quantities{parameterScanNum}(iii,1)=outputs.flux;
    quantities{parameterScanNum}(iii,2)=outputs.flow;
    quantities{parameterScanNum}(iii,3)=toc;
    quantities{parameterScanNum}(iii,4)=totalNNZ;
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
    problem = assembleMatrix(resolutionParameters, nu, E, geometryParameters, ...
        ExB_alpha_option_1,zeta_preconditioner_option_1,xi_preconditioner_option_1);
    if abs(solutionMethod)>2
        preconditioner = assembleMatrix(resolutionParameters, nu, E, geometryParameters, ...
            ExB_alpha_option_2,zeta_preconditioner_option_2,xi_preconditioner_option_2);
    else
        preconditioner=problem;
    end
    tic
    [solution, totalNNZ] = solver(problem.matrix, problem.rhs, preconditioner.matrix, solutionMethod);
    outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
    quantities{parameterScanNum}(iii,1)=outputs.flux;
    quantities{parameterScanNum}(iii,2)=outputs.flow;
    quantities{parameterScanNum}(iii,3)=toc;
    quantities{parameterScanNum}(iii,4)=totalNNZ;
end
parameterScanNum = parameterScanNum+1;

maxs=ones(numQuantities,1)*(-Inf);
mins=ones(numQuantities,1)*(Inf);
for iParameter = 1:numParameters
    maxs = max([maxs, quantities{iParameter}(:,:)'],[],2);
    mins = min([mins, quantities{iParameter}(:,:)'],[],2);
end

figure(16)
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
