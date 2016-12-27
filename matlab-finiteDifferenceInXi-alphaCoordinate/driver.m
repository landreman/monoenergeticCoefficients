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
nu = 0.1;
resolutionParameters = struct(...
    'Nalpha', 13,...
    'Nzeta', 15,...
    'Nxi', 16,...
    'includeConstraint', true);
% Should get L11 = -0.0329,  L12 = 0.011
%}
    
%{
nu = 0.01;
resolutionParameters = struct(...
    'Nalpha', 15,...
    'Nzeta', 15,...
    'Nxi', 30,...
    'includeConstraint', true);
% Should get L11 = -0.0398,  L12 = -0.238
%}

nu = 0.001;
resolutionParameters = struct(...
    'Nalpha', 27,...
    'Nzeta', 29,...
    'Nxi', 70,...
    'includeConstraint', true);
% Should get L11 = -0.09,  L12 = -1.09


%E=0.1;
E=0;
    
solutionMethod = 4;
% 1 = sparse direct solver (backslash)
% 2 = sparse direct solver (explicit LU decomposition, so memory can be monitored.)
% 3 = GMRES with no preconditioning
% 4 = GMRES with preconditioning
% -4 = same as 4 but no plot.


discretizationParameters = struct(...
    'alpha_derivative_option', 8, ...
    'zeta_derivative_option', 8, ...
    'alpha_interpolation_stencil', 4, ...
    'xi_derivative_option', 3, ...
    'pitch_angle_scattering_option', 3, ...
    'xi_quadrature_option',3);

preconditionerDiscretizationParameters = struct(...
    'alpha_derivative_option', 0, ...
    'zeta_derivative_option', 4, ...
    'alpha_interpolation_stencil', 1, ...
    'xi_derivative_option', 2, ...
    'pitch_angle_scattering_option', 2, ...
    'xi_quadrature_option',3);


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

tic
fprintf('Building main matrix.\n')
problem = assembleMatrix(resolutionParameters, nu, E, geometryParameters,discretizationParameters);
fprintf('Building preconditioner matrix.\n')
preconditioner = assembleMatrix(resolutionParameters, nu, E, geometryParameters,preconditionerDiscretizationParameters);
fprintf('Time to assemble matrix: %g sec\n',toc)

figure(1)
clf
numRows=1;
numCols=2;

subplot(numRows,numCols,1)
spy(problem.matrix)
title('Real matrix')

subplot(numRows,numCols,2)
spy(preconditioner.matrix)
title('Preconditioner matrix')


figure(2)
clf
contourf(problem.zeta2D, problem.theta2D, problem.B, 21, 'EdgeColor', 'none')
hold on
plot(problem.zeta2D, problem.theta2D, '.k')
colorbar
xlabel('\zeta')
ylabel('\theta')
title('B')

%{
fullMatrix = full(problem.matrix);
tic
cond_result = cond(fullMatrix);
fprintf('cond: %e   Time for cond: %g\n',cond_result,toc)
%}
%{
tic
condest_result = condest(problem.matrix);
fprintf('condest: %e   Time for condest: %g\n',condest_result,toc)
%}
%{
tic
rcond_result = rcond(fullMatrix);
fprintf('rcond: %e   Time for rcond: %g\n',rcond_result,toc)
%fprintf('cond=%g  condest=%g   rcond=%g\n',cond_result,condest_result,rcond_result)
%}

tic
fprintf('Beginning solve.\n')
[solution, totalNNZ] = solver(problem.matrix, problem.rhs, preconditioner.matrix, solutionMethod);
fprintf('Done. Time for solve: %g sec\n',toc)

if resolutionParameters.includeConstraint
    fprintf('Source: %g  (Should be within machine precision of 0)\n',solution(end))
end
outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
fprintf('Diagnostics: flux (L11) = %g, flow (L21) = %g.\n',outputs.flux,outputs.flow);

return


Legendres=ones(resolutionParameters.Nxi);
xi = linspace(-1,1,resolutionParameters.Nxi);
for L=1:(resolutionParameters.Nxi-1)
    temp=legendre(L,xi);
    Legendres(:,L+1) = temp(1,:)';
end
% In the 'Legendres' matrix, the first index (row) is xi, the 2nd index (col) is L

numRows=2;
numCols=3;
numContours=25;

for ialpha=1:resolutionParameters.Nalpha
    subplotNum=mod(ialpha-1,numRows*numCols)+1;
    if subplotNum==1
        figure('Position',[1,1,1600,800])
    end
    subplot(numRows,numCols,subplotNum)
    toPlot=zeros(resolutionParameters.Nxi,resolutionParameters.Nzeta);
    for izeta=1:resolutionParameters.Nzeta
        temp = solution(getIndex(ialpha,izeta,1:resolutionParameters.Nxi,resolutionParameters));
        toPlot(:,izeta)=Legendres*temp(:);
    end
    contourf(problem.zeta,xi,toPlot,numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('xi')
    title(['alpha=',num2str(problem.alpha(ialpha))])
end
