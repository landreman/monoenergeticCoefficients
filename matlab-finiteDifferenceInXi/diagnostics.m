function outputStruct = diagnostics(resolutionParameters, geometryParameters, problem, solution)

% This subroutine computes several integrals of the solution which are of interest.

Ntheta = resolutionParameters.Ntheta;
Nzeta = resolutionParameters.Nzeta;
Nxi = resolutionParameters.Nxi;

G = geometryParameters.G;
I = geometryParameters.I;

B = problem.B;
dBdtheta = problem.dBdtheta;
dBdzeta = problem.dBdzeta;
zetaWeights = problem.zetaWeights;
thetaWeights = problem.thetaWeights;
xiWeights = problem.xiWeights;
xi = problem.xi;

flux = 0;
flow = 0;

spatialPart = (1./(B.^3)) .* (G * dBdtheta - I * dBdzeta);
for itheta = 1:Ntheta
    for izeta = 1:Nzeta
        indices = getIndex(itheta,izeta,1:Nxi,resolutionParameters);

        % Compute flux:
        flux = flux + thetaWeights(itheta)*zetaWeights(izeta)*spatialPart(itheta,izeta)...
            *(xiWeights') * ((1+xi.^2) .* solution(indices));
        
        % Compute flow:
        flow = flow + thetaWeights(itheta)*zetaWeights(izeta) / B(itheta,izeta)...
            *(xiWeights') * (xi .* solution(indices));
    end
end

VPrime = thetaWeights' * (1./B.^2) * zetaWeights;
flow = flow * 2 / (sqrt(pi)*G*VPrime);
flux = -2/(sqrt(pi)*G*G*VPrime)*flux;

outputStruct = struct('flux',flux,'flow',flow);

end
