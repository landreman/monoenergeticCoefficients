clear

geometryParameters = struct(...
    'epsilon_t',-0.07053,...
    'epsilon_h',0.05067,...
    'iota',0.4542,...
    'G',3.7481,...
    'I',0,...
    'Nperiods',10,...
    'helicity_l',2,...
    'axisymmetric',false);


NFourierConverged = 20;
%NFouriers = 5:25;
NFouriers = 12:80;
%NFouriers = round(linspace(12,100,25));
%NFouriers = round(linspace(60,200,20));

NxiConverged = 50;
%Nxis = 15:25;
Nxis = round(linspace(45,90,5));

%{
% Collisionality, using the same normalization as in SFINCS:
nuPrime = 0.1;

% For comparison with Hakan's convention, scale the collisionality by nu_D(1):
Psi_Chandra = (erf(1) - 2/sqrt(pi)*exp(-1)) / 2;
nuD = 3*sqrt(pi)/4*(erf(1) - Psi_Chandra);
nu = nuPrime * nuD;
%}
nu = 0.001;

includeConstraint = true;
solutionMethod = 2;

quantitiesToRecord = {'flux','flow','time for LU','nnz(L)+nnz(U)'};

linespecs = {'.-b','.-r','.-g','.-m','.:c','.-r','.:k','.:b','.-m'};

parametersToVary = {'NFourier','N\xi'};
abscissae = {NFouriers, Nxis};
convergeds = {NFourierConverged, NxiConverged};

numQuantities = numel(quantitiesToRecord);
numParameters = numel(parametersToVary);
quantities = cell(numParameters,1);
quantities{1} = zeros(numel(NFouriers), numQuantities);
quantities{2} = zeros(numel(Nxis), numQuantities);
parameterScanNum = 1;

numRuns = numel(NFouriers) + numel(Nxis);
runNum = 1;

% Vary Ntheta, keeping other numerical parameters fixed.
for iii = 1:numel(NFouriers)
    resolutionParameters = struct(...
        'NFourier', NFouriers(iii),...
        'NFourier2', NFouriers(iii)*2-1,...
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
        'NFourier', NFourierConverged,...
        'NFourier2', NFourierConverged*2-1,...
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
