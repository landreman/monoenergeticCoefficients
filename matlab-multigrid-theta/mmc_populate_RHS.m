function mmc_populate_RHS()

global rhs G I levels
global Ntheta Nzeta Nxi

resolutionParameters = struct('Ntheta',Ntheta,'Nzeta',Nzeta,'Nxi',Nxi,'matrixSize',levels(1).matrixSize);

rhs = zeros(levels(1).matrixSize,1);
xi = levels(1).xi;
spatialPart = (1./levels(1).B) .* (G * levels(1).dBdtheta - I * levels(1).dBdzeta);
for itheta=1:Ntheta
    for izeta = 1:Nzeta
        indices = mmc_get_index(itheta,izeta,1:Nxi,resolutionParameters);
        rhs(indices) = (1 + xi.^2) * spatialPart(itheta,izeta);
    end
end


end
