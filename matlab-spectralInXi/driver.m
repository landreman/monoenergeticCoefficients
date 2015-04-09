clear

geometryParameters = struct(...
    'epsilon_t',-0.07053,...
    'epsilon_h',0.05067,...
    'iota',0.4542,...
    'G',3.7481,...
    'I',0,...
    'Nperiods',10,...
    'helicity_l',2);

resolutionParameters = struct(...
    'Ntheta', 13,...
    'Nzeta', 15,...
    'Nxi', 16);

% Collisionality:
% (Typically nu << 1, so advection dominates diffusion.)
nu = 0.1;

solutionMethod = 2;
% 1 = sparse direct solver (backslash)
% 2 = GMRES with no preconditioning

% **********************************************************
% End of input parameters.
% **********************************************************

tic
problem = assembleMatrix(resolutionParameters, nu, geometryParameters);
fprintf('Time to assemble matrix: %g sec\n',toc)

figure(1)
clf
spy(problem.matrix)

figure(2)
clf
contourf(problem.zeta2D, problem.theta2D, problem.B, 21, 'EdgeColor', 'none')
hold on
plot(problem.zeta2D, problem.theta2D, '.k')
colorbar
xlabel('\zeta')
ylabel('\theta')
title('B')

tic
fprintf('Beginning solve.\n')
solution = solver(problem.matrix, problem.rhs, solutionMethod);
fprintf('Done. Time for solve: %g sec\n',toc)

fprintf('Source: %g  (Should be within machine precision of 0)\n',solution(end))
outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
fprintf('Diagnostics: flux (L11) = %g, flow (L21) = %g.\n',outputs.flux,outputs.flow);
