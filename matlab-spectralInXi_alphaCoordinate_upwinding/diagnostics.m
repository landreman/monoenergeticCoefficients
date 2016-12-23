function outputStruct = diagnostics(resolutionParameters, geometryParameters, problem, solution)

% This subroutine computes several integrals of the solution which are of interest.

Nalpha = resolutionParameters.Nalpha;
Nzeta = resolutionParameters.Nzeta;

G = geometryParameters.G;
I = geometryParameters.I;

B = problem.B;
dBdtheta = problem.dBdtheta;
dBdzeta = problem.dBdzeta;
zetaWeights = problem.zetaWeights;
alphaWeights = problem.alphaWeights;

flux = 0;
flow = 0;

spatialPart = (1./(B.^3)) .* (G * dBdtheta - I * dBdzeta);
for ialpha = 1:Nalpha
    % Compute flux:
    L=0;
    indices = getIndex(ialpha,1:Nzeta,L+1,resolutionParameters);
    flux = flux + (8/3)*alphaWeights(ialpha)*(zetaWeights' .* spatialPart(ialpha,:))*solution(indices);
    
    L=2;
    indices = getIndex(ialpha,1:Nzeta,L+1,resolutionParameters);
    flux = flux + (4/15)*alphaWeights(ialpha)*(zetaWeights' .* spatialPart(ialpha,:))*solution(indices);

    % Compute flow:
    L=1;
    indices = getIndex(ialpha,1:Nzeta,L+1,resolutionParameters);
    flow = flow + alphaWeights(ialpha)*(zetaWeights' ./ B(ialpha,:))*solution(indices);

end

VPrime = alphaWeights' * (1./B.^2) * zetaWeights;
flow = flow * 4 / (3*sqrt(pi)*G*VPrime);
flux = -2/(sqrt(pi)*G*G*VPrime)*flux;

outputStruct = struct('flux',flux,'flow',flow);

end
