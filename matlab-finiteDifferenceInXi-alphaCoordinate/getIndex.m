function index = getIndex(ialpha,izeta,ixi,resolutionParameters)

Nalpha = resolutionParameters.Nalpha;
Nzeta = resolutionParameters.Nzeta;
Nxi = resolutionParameters.Nxi;

% Validation:
assert(all(ialpha>=1))
assert(all(izeta>=1))
assert(all(ixi>=1))
assert(all(ialpha <= Nalpha))
assert(all(izeta <= Nzeta))
assert(all(ixi <= Nxi))

% The key formula:
index = (ixi-1)*Nalpha*Nzeta ...
    + (ialpha-1)*Nzeta ...
    + izeta;

end