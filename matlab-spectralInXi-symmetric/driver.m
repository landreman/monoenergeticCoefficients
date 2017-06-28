%clear

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
nu = 1e-1;
E=0;
% Results should be flux (L11) = -0.0329, flow (L21) = 0.0111
    
%{
resolutionParameters = struct(...
    'Ntheta', 15,...
    'Nzeta', 15,...
    'Nxi', 30);
nu = 1e-2;
E=0;
% Results should be flux (L11) = -0.0392731, flow (L21) = -0.237442.
%}

% Collisionality:
% (Typically nu << 1, so advection dominates diffusion.)

%{
% Collisionality, using the same normalization as in SFINCS:
nuPrime = 0.1;

% For comparison with Hakan's convention, scale the collisionality by nu_D(1):
Psi_Chandra = (erf(1) - 2/sqrt(pi)*exp(-1)) / 2;
nuD = 3*sqrt(pi)/4*(erf(1) - Psi_Chandra);
nu = nuPrime * nuD;
fprintf('For SFINCS nuPrime=%g, the equivalent nu for the monoenergetic code=%g\n',nuPrime,nu)
%}



solutionMethod = 2;
% 1 = sparse direct solver (backslash)
% 2 = sparse direct solver (explicit LU decomposition, so memory can be monitored.)
% 3 = GMRES with no preconditioning

% **********************************************************
% End of input parameters.
% **********************************************************

tic
problem = assembleMatrix(resolutionParameters, nu, E, geometryParameters);
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
[solution, totalNNZ] = solver(problem.matrix, problem.rhs, solutionMethod);
fprintf('Done. Time for solve: %g sec\n',toc)

fprintf('Sources: %g, %g  (Should be within machine precision of 0)\n',solution(end-1), solution(end))
outputs = diagnostics(resolutionParameters, geometryParameters, problem, solution);
fprintf('Diagnostics: flux (L11) = %g, flow (L21) = %g.\n',outputs.flux,outputs.flow);
