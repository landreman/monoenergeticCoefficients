function returnStruct = assembleMatrix(resolutionParameters, nu, geometryParameters)

NFourier = resolutionParameters.NFourier;
Nxi = resolutionParameters.Nxi;
NFourier2 = resolutionParameters.NFourier2;

matrixSize = NFourier2*Nxi;
if resolutionParameters.includeConstraint
    matrixSize = matrixSize + 1;
end

% Choose (m,n) pairs to use:
[ms,ns,Binv_vec] = chooseFourierModes(geometryParameters,resolutionParameters);

[ddtheta,ddzeta] = buildFourierDifferentiationMatrices(ms,ns);


% Make a vector in mn-space that represents |B|:
Bvec = zeros(NFourier2,1);
% Constant term:
Bvec(1) = 1;
% Toroidal mode:
mask = (ms==1 & ns==0);
index = find(mask);
if numel(index) ~= 1
    error('Error! Did not find the toroidal mode!\n')
end
Bvec(index) = geometryParameters.epsilon_t;
% Helical mode:
mask = (ms==geometryParameters.helicity_l & ns==geometryParameters.Nperiods);
index = find(mask);
if numel(index) ~= 1
    error('Error! Did not find the helical mode!\n')
end
Bvec(index) = geometryParameters.epsilon_h;

BMatrix = buildFourierConvolutionMatrix(ms,ns,Bvec);
BinvMatrix = buildFourierConvolutionMatrix(ms,ns,Binv_vec);
iota = geometryParameters.iota;

estimated_nnz = ...
    NFourier2*NFourier2*Nxi*2 ... % streaming term (dense in Fourier, 2 off-diagonals in xi)
    + NFourier2*NFourier2*Nxi*2 ... % ddxi term (dense in Fourier, 2 off-diagonals in xi)
    + NFourier2*(Nxi-1) ...   % collision term (diagonal, 0 when L=0)
    + 1 ... % Source
    + NFourier2;  % Constraint

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
% Add streaming terms:
% -----------------------------------------
iota_ddtheta_plus_ddzeta = iota*ddtheta+ddzeta;

spatialMatrix = BMatrix*iota_ddtheta_plus_ddzeta;

for L=0:(Nxi-1)
    rowIndices = getIndex(1:NFourier2, L+1, resolutionParameters);
    
    % Super-diagonal term
    if (L<Nxi-1)
        ell = L + 1;
        colIndices = getIndex(1:NFourier2,ell+1,resolutionParameters);
        addSparseBlock(rowIndices, colIndices, (L+1)/(2*L+3)*spatialMatrix)
    end
    
    % Sub-diagonal term
    if (L>0)
        ell = L - 1;
        colIndices = getIndex(1:NFourier2,ell+1,resolutionParameters);
        addSparseBlock(rowIndices, colIndices, L/(2*L-1)*spatialMatrix)
    end
    
end

% -----------------------------------------
% Add d/dxi terms:
% -----------------------------------------


%spatialPartOfTerm = -(iota*dBdtheta(itheta,:)+dBdzeta(itheta,:))/2;
spatialPartOfTerm_vec = -(1/2)*iota_ddtheta_plus_ddzeta*Bvec;
spatialPartOfTerm = buildFourierConvolutionMatrix(ms,ns,spatialPartOfTerm_vec);
for L=0:(Nxi-1)
    rowIndices = getIndex(1:NFourier2, L+1, resolutionParameters);
    
    % Super-diagonal term
    if (L<Nxi-1)
        ell = L + 1;
        colIndices = getIndex(1:NFourier2,ell+1,resolutionParameters);
        addSparseBlock(rowIndices, colIndices, (L+1)*(L+2)/(2*L+3)*spatialPartOfTerm)
    end
    
    % Sub-diagonal term
    if (L>0)
        ell = L - 1;
        colIndices = getIndex(1:NFourier2,ell+1,resolutionParameters);
        addSparseBlock(rowIndices, colIndices, (-L)*(L-1)/(2*L-1)*spatialPartOfTerm)
    end
    
end

% -----------------------------------------
% Add the diffusion (collision) term:
% -----------------------------------------

L = (1:Nxi)-1;
for imn=1:NFourier2
    indices = getIndex(imn,L+1,resolutionParameters);
    addToSparse(indices, indices, nu/2*L.*(L+1))
end

% -----------------------------------------
% Add the extra constraint:
% -----------------------------------------

if resolutionParameters.includeConstraint
    rowIndex = matrixSize;
    L = 0;
    BinvMatrixSquared = BinvMatrix*BinvMatrix;
    colIndices = getIndex(1:NFourier2,L+1,resolutionParameters);
    %addSparseBlock(rowIndex, colIndices, thetaWeights(itheta) * (zetaWeights') ./ (B(itheta,:) .^ 2))
    addSparseBlock(rowIndex, colIndices, BinvMatrixSquared(1,:))
end

% -----------------------------------------
% Add the "source" lambda:
% -----------------------------------------

if resolutionParameters.includeConstraint
    colIndex = matrixSize;
    L = 0;
    rowIndex = getIndex(1,L+1,resolutionParameters);
    addSparseBlock(rowIndex, colIndex, 1)
end


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

dBdtheta = ddtheta * Bvec;
dBdzeta = ddzeta  * Bvec;
rhs = zeros(matrixSize,1);
%spatialPart = (1./B) .* (geometryParameters.G * dBdtheta - geometryParameters.I * dBdzeta);
spatialPart = BinvMatrix * (geometryParameters.G * dBdtheta - geometryParameters.I * dBdzeta);

L=0;
indices = getIndex(1:NFourier2,L+1,resolutionParameters);
rhs(indices) = spatialPart * (4/3);

L=2;
indices = getIndex(1:NFourier2,L+1,resolutionParameters);
rhs(indices) = spatialPart * (2/3);


% -----------------------------------------
% Done. Return a structure with everything important:
% -----------------------------------------

returnStruct = struct(...
    'matrixSize',matrixSize,...
    'ddtheta',ddtheta,...
    'ddzeta',ddzeta,...
    'BinvMatrix',BinvMatrix,...
    'dBdtheta',dBdtheta,...
    'dBdzeta',dBdzeta,...
    'matrix',matrix,...
    'ms',ms,...
    'ns',ns,...
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
