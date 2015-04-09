function returnStruct = assembleMatrix(resolutionParameters, nu, geometryParameters)

Ntheta = resolutionParameters.Ntheta;
Nzeta = resolutionParameters.Nzeta;
Nxi = resolutionParameters.Nxi;

matrixSize = Ntheta*Nzeta*Nxi;

% Build theta grid:

%scheme = 0; % Uniform periodic 2nd order FD                    
scheme = 10; % Uniform periodic 4th order FD
%[theta, ddtheta, thetaWeights] = setupGrid(Ntheta,0,2*pi);
[theta, thetaWeights, ddtheta, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Ntheta, 0, 2*pi, scheme);

% Build zeta grid:
if Ntheta==1
    zeta=0;
    ddzeta=0;
    zetaWeights=2*pi;
else
    zetaMax = 2*pi/geometryParameters.Nperiods;
    %[zeta, ddzeta, zetaWeights] = setupGrid(Nzeta,0,zetaMax);
    %zetaWeights = zetaWeights * geometryParameters.Nperiods;

    %scheme = 0; % Uniform periodic 2nd order FD                    
    scheme = 10; % Uniform periodic 4th order FD
    [zeta, zetaWeights, ddzeta, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Nzeta, 0, zetaMax, scheme);
end

% Switch to cell-centered:
theta = theta + (theta(2)-theta(1))/2;
zeta = zeta + (zeta(2)-zeta(1))/2;

[zeta2D, theta2D] = meshgrid(zeta,theta);

% Build xi grid:
xi = linspace(-1,1,Nxi+1);
xi = xi(:);
xi(end)=[];
dxi = xi(2)-xi(1);
xi = xi + dxi/2; % Make cell-centered

scheme = 2; % First and last rows use one-sided 3-point stencil
%scheme = 3; % First and last rows use one-sided 2-point stencil
[~, ~, ddxi_mirror, ~] = m20121125_04_DifferentiationMatricesForUniformGrid(Nxi, min(xi), max(xi), scheme);

xiWeights = zeros(size(xi))+2/Nxi;

% collisionOperator = (d/dxi) (1-xi^2) (d/dxi)
dxi2 = dxi*dxi;
collisionOperator = ...
    diag((1-(xi(1:end-1)+dxi/2).^2)/dxi2, 1) ...
    - diag((1-(xi+dxi/2).^2)/dxi2 + (1-(xi-dxi/2).^2)/dxi2, 0) ...
    + diag((1-(xi(2:end)-dxi/2).^2)/dxi2, -1);


% Initialize the arrays for the magnetic field strength B:
iota = geometryParameters.iota;
[B, dBdtheta, dBdzeta] = geometry(theta2D, zeta2D, geometryParameters);

estimated_nnz = ...
    Ntheta*Nzeta*Nxi*4*2 ... % ddtheta term (4 off-diagonals in theta, 2 off-diagonals in xi)
    + Ntheta*Nzeta*Nxi*4*2 ...  % ddzeta term (4 off-diagonals in zeta, 2 off-diagonals in xi)
    + Ntheta*Nzeta*Nxi*2 ... % ddxi term (2 off-diagonals in xi)
    + Ntheta*Nzeta*(Nxi-1) ...   % collision term (diagonal, 0 when L=0)
    + Ntheta*Nzeta ... % constraint term
    + Ntheta*Nzeta; % source term

% Initialize some arrays that are used for building the sparse matrix:
sparseCreatorIndex=1;
sparseCreator_i=0;
sparseCreator_j=0;
sparseCreator_s=0;
resetSparseCreator()

% ***************************************************************************
% ***************************************************************************
% Now populate the matrix.
% ***************************************************************************
% ***************************************************************************

% -----------------------------------------
% Add d/dtheta terms:
% -----------------------------------------

for izeta=1:Nzeta
    thetaPartOfTerm = iota*diag(B(:,izeta))*ddtheta;
    for ixi = 1:Nxi
        indices = getIndex(1:Ntheta,izeta,ixi,resolutionParameters);
        addSparseBlock(indices, indices, xi(ixi)*thetaPartOfTerm)
    end
end

% -----------------------------------------
% Add d/dzeta terms:
% -----------------------------------------

for itheta=1:Ntheta
    zetaPartOfTerm = diag(B(itheta,:))*ddzeta;
    for ixi = 1:Nxi
        indices = getIndex(itheta, 1:Nzeta, ixi, resolutionParameters);
        addSparseBlock(indices, indices, xi(ixi)*zetaPartOfTerm)
    end
end

% -----------------------------------------
% Add mirror term:
% -----------------------------------------

xiPart = diag(1-xi.^2)*ddxi_mirror;
spatialPart = -(iota*dBdtheta+dBdzeta)/2;
for itheta=1:Ntheta
    for izeta = 1:Nzeta
        indices = getIndex(itheta, izeta, 1:Nxi, resolutionParameters);
        addSparseBlock(indices, indices, spatialPart(itheta,izeta)*xiPart)
    end
end

% -----------------------------------------
% Add the diffusion (collision) term:
% -----------------------------------------

for itheta=1:Ntheta
    for izeta=1:Nzeta
        indices = getIndex(itheta,izeta,1:Nxi,resolutionParameters);
        addSparseBlock(indices, indices, -nu/2*collisionOperator)
    end
end

%{
% -----------------------------------------
% Add the extra constraint:
% -----------------------------------------

rowIndex = matrixSize;
L = 0;
for itheta=1:Ntheta
    colIndices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    addSparseBlock(rowIndex, colIndices, thetaWeights(itheta) * (zetaWeights') ./ (B(itheta,:) .^ 2))
end

% -----------------------------------------
% Add the "source" lambda:
% -----------------------------------------

colIndex = matrixSize;
L = 0;
for itheta=1:Ntheta
    rowIndices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    addSparseBlock(rowIndices, colIndex, ones(Nzeta,1))
end
%}

% -----------------------------------------
% Finalize the matrix
% -----------------------------------------
matrix = createSparse();

fprintf('matrixSize: %g,  nnz: %g,  estimated nnz: %g,  sparsity: %g\n',...
    matrixSize,nnz(matrix),estimated_nnz, nnz(matrix)/(matrixSize^2))

% ***************************************************************************
% ***************************************************************************
% Done with assembling the matrix.
% ***************************************************************************
% ***************************************************************************

% -----------------------------------------
% Build the right-hand side
% -----------------------------------------

rhs = zeros(matrixSize,1);
spatialPart = (1./B) .* (geometryParameters.G * dBdtheta - geometryParameters.I * dBdzeta);
for itheta = 1:Ntheta
    for izeta = 1:Nzeta
        indices = getIndex(itheta,izeta,1:Nxi,resolutionParameters);
        rhs(indices) = (1 + xi.^2) * spatialPart(itheta,izeta);
    end
end


% -----------------------------------------
% Done. Return a structure with everything important:
% -----------------------------------------

returnStruct = struct(...
    'matrixSize',matrixSize,...
    'theta',theta,...
    'zeta',zeta,...
    'xi',xi,...
    'theta2D',theta2D,...
    'zeta2D',zeta2D,...
    'thetaWeights',thetaWeights,...
    'zetaWeights',zetaWeights,...
    'xiWeights',xiWeights,...
    'B',B,...
    'dBdtheta',dBdtheta,...
    'dBdzeta',dBdzeta,...
    'matrix',matrix,...
    'rhs',rhs...
    );
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions below are all utilities for building sparse matrices:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function resetSparseCreator()
        sparseCreatorIndex=1;
        sparseCreator_i=zeros(estimated_nnz,1);
        sparseCreator_j=zeros(estimated_nnz,1);
        sparseCreator_s=zeros(estimated_nnz,1);
    end

    function addToSparse(i,j,s)
        % Adds values to the sparse matrix.
        n=numel(i);
        if n ~= numel(j)
            error('Error A');
        end
        if n ~= numel(s)
            error('Error B');
        end
        if any(i<1)
            error('Error Q: i<1');
        end
        if any(j<1)
            error('Error Q: j<1');
        end
        sparseCreator_i(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = i;
        sparseCreator_j(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = j;
        sparseCreator_s(sparseCreatorIndex:(sparseCreatorIndex+n-1)) = s;
        sparseCreatorIndex = sparseCreatorIndex+n;
        if sparseCreatorIndex > estimated_nnz
            fprintf('Warning! estimated_nnz is too small. Increase predictedNNZ.\n')
        end
    end

    function addSparseBlock(rowIndices, colIndices, block)
        % Adds a block to the sparse matrix.
        % rowIndices and colIndices should be vectors.
        % numel(rowIndices) should equal the number of rows in 'block'.
        % numel(colIndices) should equal the number of columns in 'block'.
        s=size(block);
        if (s(1) ~= numel(rowIndices)) || (s(2) ~= numel(colIndices))
            error('Error in addSparseBlock!  Input sizes are not consistent.')
        end
        [rows, cols, values] = find(block);
        addToSparse(rowIndices(rows),colIndices(cols),values)
    end

    function sparseMatrix = createSparse()
        % After you are done adding elements to the sparse matrix using
        % addToSparse() and addSparseBlock(), call this function to
        % finalize the matrix.
        sparseMatrix = sparse(sparseCreator_i(1:(sparseCreatorIndex-1)), sparseCreator_j(1:(sparseCreatorIndex-1)), sparseCreator_s(1:(sparseCreatorIndex-1)), matrixSize, matrixSize);
        resetSparseCreator()
    end

end
