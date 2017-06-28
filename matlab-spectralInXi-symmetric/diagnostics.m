function outputStruct = diagnostics(resolutionParameters, geometryParameters, problem, solution)

% This subroutine computes several integrals of the solution which are of interest.

Ntheta = resolutionParameters.Ntheta;
Nzeta = resolutionParameters.Nzeta;

G = geometryParameters.G;
I = geometryParameters.I;

B = problem.B;
dBdtheta = problem.dBdtheta;
dBdzeta = problem.dBdzeta;
zetaWeights = problem.zetaWeights;
thetaWeights = problem.thetaWeights;
L_scaling = problem.L_scaling;

flux = 0;
flow = 0;

spatialPart = (1./(B.^3)) .* (G * dBdtheta - I * dBdzeta);
for itheta = 1:Ntheta
    % Compute flux:
    L=0;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    flux = flux + L_scaling(L+1)*(8/3)*thetaWeights(itheta)*(zetaWeights' .* spatialPart(itheta,:))*solution(indices);
    
    L=2;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    flux = flux + L_scaling(L+1)*(4/15)*thetaWeights(itheta)*(zetaWeights' .* spatialPart(itheta,:))*solution(indices);

    % Compute flow:
    L=1;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    flow = flow + L_scaling(L+1)*thetaWeights(itheta)*(zetaWeights' ./ B(itheta,:))*solution(indices);

end

%VPrime had been int int 1/B^2. Now it's int int sqrt_g.
VPrime_old = thetaWeights' * (1./(B.*B)) * zetaWeights;
flow = flow * 4 / (3*sqrt(pi)*G*VPrime_old);
flux = -2/(sqrt(pi)*G*G*VPrime_old)*flux;

outputStruct = struct('flux',flux,'flow',flow);

end
