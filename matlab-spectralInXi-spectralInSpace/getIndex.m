function index = getIndex(imn,ixi,resolutionParameters)

NFourier2 = resolutionParameters.NFourier2;
Nxi = resolutionParameters.Nxi;

% Validation:
assert(all(imn>=1))
assert(all(ixi>=1))
assert(all(imn <= NFourier2))
assert(all(ixi <= Nxi))

% The key formula:
index = (ixi-1)*NFourier2 ...
    + imn;

end