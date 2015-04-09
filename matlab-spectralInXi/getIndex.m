function index = getIndex(itheta,izeta,ixi,resolutionParameters)

Ntheta = resolutionParameters.Ntheta;
Nzeta = resolutionParameters.Nzeta;
Nxi = resolutionParameters.Nxi;

% Validation:
assert(all(itheta>=1))
assert(all(izeta>=1))
assert(all(ixi>=1))
assert(all(itheta <= Ntheta))
assert(all(izeta <= Nzeta))
assert(all(ixi <= Nxi))

% The key formula:
index = (ixi-1)*Ntheta*Nzeta ...
    + (itheta-1)*Nzeta ...
    + izeta;

end