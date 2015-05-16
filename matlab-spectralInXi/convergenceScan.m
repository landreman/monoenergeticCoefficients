clear

geometryParameters = struct(...
    'epsilon_t',-0.07053,...
    'epsilon_h',0.05067,...
    'iota',0.4542,...
    'G',3.7481,...
    'I',0,...
    'Nperiods',10,...
    'helicity_l',2);


NthetaConverged = 13;
Nthetas = 9:2:29;

NzetaConverged = 13;
Nzetas = 9:2:29;

NxiConverged = 16;
Nxis = 15:25;

% Collisionality, using the same normalization as in SFINCS:
nuPrime = 0.1;

% For comparison with Hakan's convention, scale the collisionality by nu_D(1):
Psi_Chandra = (erf(1) - 2/sqrt(pi)*exp(-1)) / 2;
nuD = 3*sqrt(pi)/4*(erf(1) - Psi_Chandra);
nu = nuPrime * nuD;

includeConstraint = true;
solutionMethod = 2;

quantitiesToRecord = {'flux','flow','time for LU','nnz(L)+nnz(U)'};

linespecs = {'.-b','.-r','.-g','.-m','.:c','.-r','.:k','.:b','.-m'};

parametersToVary = {'N\theta','N\zeta','N\xi'};
abscissae = {Nthetas, Nzetas, Nxis};
convergeds = {NthetaConverged, NzetaConverged, NxiConverged};

numQuantities = numel(quantitiesToRecord);
numParameters = numel(parametersToVary);
quantities = cell(numParameters,1);
quantities{1} = zeros(numel(Nthetas), numQuantities);
quantities{2} = zeros(numel(Nzetas), numQuantities);
quantities{3} = zeros(numel(Nxis), numQuantities);
parameterScanNum = 1;

numRuns = numel(Nthetas)+numel(Nzetas)+numel(Nxis);
runNum = 1;

% Vary Ntheta, keeping other numerical parameters fixed.
for iii = 1:numel(Nthetas)
    resolutionParameters = struct(...
        'Ntheta', Nthetas(iii),...
        'Nzeta', NzetaConverged,...
        'Nxi', NxiConverged,...
        'includeConstraint', includeConstraint);
    fprintf('Beginning solve %d of %d.\n',runNum,numRuns); runNum = runNum+1;
    problem = assembleMatrix(resolutionParameters, nu, geometryParameters);
    tic
    [solution, totalNNZ] = solver(problem.matrix, problem.rhs, solutionMethod);
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
        'Ntheta', NthetaConverged,...
        'Nzeta', Nzetas(iii),...
        'Nxi', NxiConverged,...
        'includeConstraint', includeConstraint);
    fprintf('Beginning solve %d of %d.\n',runNum,numRuns); runNum = runNum+1;
    problem = assembleMatrix(resolutionParameters, nu, geometryParameters);
    tic
    [solution, totalNNZ] = solver(problem.matrix, problem.rhs, solutionMethod);
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
        'Ntheta', NthetaConverged,...
        'Nzeta', NzetaConverged,...
        'Nxi', Nxis(iii),...
        'includeConstraint', includeConstraint);
    fprintf('Beginning solve %d of %d.\n',runNum,numRuns); runNum = runNum+1;
    problem = assembleMatrix(resolutionParameters, nu, geometryParameters);
    tic
    [solution, totalNNZ] = solver(problem.matrix, problem.rhs, solutionMethod);
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

figure(3)
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
