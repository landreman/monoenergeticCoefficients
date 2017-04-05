function mmc_diagnostics()

% This subroutine computes several integrals of the solution which are of interest.

global Ntheta Nzeta Nxi
global levels solution
global flux flow VPrime G I

B = levels(1).B;
dBdtheta = levels(1).dBdtheta;
dBdzeta = levels(1).dBdzeta;
zetaWeights = levels(1).zetaWeights;
thetaWeights = levels(1).thetaWeights;
xiWeights = levels(1).xiWeights;
xi = levels(1).xi;
resolutionParameters = struct('Ntheta',Ntheta,'Nzeta',Nzeta,'Nxi',Nxi,'matrixSize',levels(1).matrixSize);

flux = 0;
flow = 0;

spatialPart = (1./(B.^3)) .* (G * dBdtheta - I * dBdzeta);
for itheta = 1:Ntheta
    for izeta = 1:Nzeta
        indices = mmc_get_index(itheta,izeta,1:Nxi,resolutionParameters);
        
        % Compute flux:
        flux = flux + thetaWeights(itheta) * zetaWeights(izeta) * spatialPart(itheta,izeta) ...
            * (xiWeights')*((1+xi.^2) .* solution(indices));
        
        % Compute flow:
        flow = flow + thetaWeights(itheta) * zetaWeights(izeta) / B(itheta,izeta) ...
            * (xiWeights') * (xi .* solution(indices));
    end
end

VPrime = thetaWeights' * (1./B.^2) * zetaWeights;
flow = flow * 2 / (sqrt(pi)*G*VPrime);
flux = -2/(sqrt(pi)*G*G*VPrime)*flux;

end
