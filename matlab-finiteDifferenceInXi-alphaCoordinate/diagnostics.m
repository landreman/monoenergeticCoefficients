function outputStruct = diagnostics(resolutionParameters, geometryParameters, problem, solution)

% This subroutine computes several integrals of the solution which are of interest.

Nalpha = resolutionParameters.Nalpha;
Nzeta = resolutionParameters.Nzeta;
Nxi = resolutionParameters.Nxi;

G = geometryParameters.G;
I = geometryParameters.I;

B = problem.B;
dBdtheta = problem.dBdtheta;
dBdzeta = problem.dBdzeta;
zetaWeights = problem.zetaWeights;
alphaWeights = problem.alphaWeights;
xiWeights = problem.xiWeights;
xi = problem.xi;

flux = 0;
flow = 0;

spatialPart = (1./(B.^3)) .* (G * dBdtheta - I * dBdzeta);
for ialpha = 1:Nalpha
    for izeta = 1:Nzeta
        indices = getIndex(ialpha,izeta,1:Nxi,resolutionParameters);
        
        % Compute flux:
        flux = flux + alphaWeights(ialpha) * zetaWeights(izeta) * spatialPart(ialpha,izeta) ...
            * (xiWeights')*((1+xi.^2) .* solution(indices));
        
        % Compute flow:
        flow = flow + alphaWeights(ialpha) * zetaWeights(izeta) / B(ialpha,izeta) ...
            * (xiWeights') * (xi .* solution(indices));
    end
end

VPrime = alphaWeights' * (1./B.^2) * zetaWeights;
flow = flow * 2 / (sqrt(pi)*G*VPrime);
flux = -2/(sqrt(pi)*G*G*VPrime)*flux;

outputStruct = struct('flux',flux,'flow',flow);

end
