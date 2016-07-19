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

NFourier = 10;
Nxi = 16;

resolutionParameters = struct(...
    'NFourier', NFourier,...
    'NFourier2', NFourier*2-1,...
    'Nxi', Nxi,...
    'includeConstraint', true);

% Collisionality:
% (Typically nu << 1, so advection dominates diffusion.)
nu = 0.1;

solutionMethod = 2;
% 1 = sparse direct solver (backslash)
% 2 = sparse direct solver (explicit LU decomposition, so memory can be monitored.)
% 3 = GMRES with no preconditioning

% **********************************************************
% End of input parameters.
% **********************************************************

tic
problem = assembleMatrix(resolutionParameters, nu, geometryParameters);
fprintf('Time to assemble matrix: %g sec\n',toc)

figure(1)
clf
spy(problem.matrix)

%{
figure(2)
clf
contourf(problem.zeta2D, problem.theta2D, problem.B, 21, 'EdgeColor', 'none')
hold on
plot(problem.zeta2D, problem.theta2D, '.k')
colorbar
xlabel('\zeta')
ylabel('\theta')
title('B')
%}

tic
fprintf('Beginning solve.\n')
[solution, totalNNZ] = solver(problem.matrix, problem.rhs, solutionMethod);
fprintf('Done. Time for solve: %g sec\n',toc)

if resolutionParameters.includeConstraint
    fprintf('Source: %g  (Should be within machine precision of 0)\n',solution(end))
end
outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
fprintf('Diagnostics: flux (L11) = %g, flow (L21) = %g.\n',outputs.flux,outputs.flow);



mmax = max(problem.ms);
nmax = max(abs(problem.ns))/geometryParameters.Nperiods;

numPlots = min(Nxi,16);
numCols = ceil(sqrt(numPlots));
numRows = ceil(numPlots/numCols);

figure(4)
clf
for whichPlot = 1:numPlots
    subplot(numRows,numCols,whichPlot)
    L=whichPlot-1;
    data = NaN*zeros(mmax+1,2*nmax+1);
    f = solution(getIndex(1:resolutionParameters.NFourier2, L+1, resolutionParameters));
    data(1,nmax+1) = abs(f(1));
    for imn=2:resolutionParameters.NFourier
        data(problem.ms(imn)+1, problem.ns(imn)/geometryParameters.Nperiods+nmax+1) = sqrt(f(imn)^2 + f(imn+NFourier-1)^2);
    end
    imagesc([-nmax,nmax],[0,mmax],log10(data))
    colorbar
    xlabel('n')
    ylabel('m')
    title(['L=',num2str(L)])
end