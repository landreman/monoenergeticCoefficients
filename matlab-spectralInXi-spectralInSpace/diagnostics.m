function outputStruct = diagnostics(resolutionParameters, geometryParameters, problem, solution)

% This subroutine computes several integrals of the solution which are of interest.

NFourier2 = resolutionParameters.NFourier2;

G = geometryParameters.G;
I = geometryParameters.I;

ms = problem.ms;
ns = problem.ns;
BinvMatrix = problem.BinvMatrix;
dBdtheta = problem.dBdtheta; % This is a vector.
dBdzeta = problem.dBdzeta; % This is a vector.
BinvMatrixSquared = BinvMatrix*BinvMatrix;

% Compute flux:
flux = 0;
%spatialPart = (1./(B.^3)) .* (G * dBdtheta - I * dBdzeta);
spatialPart = buildFourierConvolutionMatrix(ms, ns, BinvMatrixSquared*BinvMatrix * (G * dBdtheta - I * dBdzeta));

L=0;
indices = getIndex(1:NFourier2,L+1,resolutionParameters);
tempVec = spatialPart * solution(indices);
flux = flux + (8/3)*4*pi*pi*tempVec(1);

L=2;
indices = getIndex(1:NFourier2,L+1,resolutionParameters);
tempVec = spatialPart * solution(indices);
flux = flux + (4/15)*4*pi*pi*tempVec(1);

% Compute flow:
L=1;
indices = getIndex(1:NFourier2,L+1,resolutionParameters);
tempVec = BinvMatrix * solution(indices);
flow = 4*pi*pi*tempVec(1);



%VPrime = thetaWeights' * (1./B.^2) * zetaWeights;
VPrime = BinvMatrixSquared(1,1)*4*pi*pi;
fprintf('VPrime: %g\n',VPrime)
flow = flow * 4 / (3*sqrt(pi)*G*VPrime);
flux = -2/(sqrt(pi)*G*G*VPrime)*flux;

outputStruct = struct('flux',flux,'flow',flow);

end
