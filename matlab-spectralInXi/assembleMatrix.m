function returnStruct = assembleMatrix(resolutionParameters, nu, geometryParameters)

Ntheta = resolutionParameters.Ntheta;
Nzeta = resolutionParameters.Nzeta;
Nxi = resolutionParameters.Nxi;

matrixSize = Ntheta*Nzeta*Nxi;
if resolutionParameters.includeConstraint
    matrixSize = matrixSize + 1;
end

% Build theta grid:
[theta, ddtheta, thetaWeights] = setupGrid(Ntheta,0,2*pi);

% Build zeta grid:
if Nzeta==1
    zeta=0;
    ddzeta=0;
    zetaWeights=2*pi;
else
    zetaMax = 2*pi/geometryParameters.Nperiods;
    [zeta, ddzeta, zetaWeights] = setupGrid(Nzeta,0,zetaMax);
    zetaWeights = zetaWeights * geometryParameters.Nperiods;
end

[zeta2D, theta2D] = meshgrid(zeta,theta);

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

Ls = (1:Nxi)-1;
%rowScaling = (2*Ls+1) .^ (1.5);
rowScaling = ones(size(Ls));

% -----------------------------------------
% Add d/dtheta terms:
% -----------------------------------------

for izeta=1:Nzeta
    thetaPartOfTerm = iota*diag(B(:,izeta))*ddtheta;
    for L=0:(Nxi-1)
        rowIndices = getIndex(1:Ntheta,izeta,L+1,resolutionParameters);
        
        % Super-diagonal term
        if (L<Nxi-1)
            ell = L + 1;
            colIndices = getIndex(1:Ntheta,izeta,ell+1,resolutionParameters);
            addSparseBlock(rowIndices, colIndices, (L+1)/(2*L+3)*thetaPartOfTerm*rowScaling(L+1))
        end
        
        % Sub-diagonal term
        if (L>0)
            ell = L - 1;
            colIndices = getIndex(1:Ntheta,izeta,ell+1,resolutionParameters);
            addSparseBlock(rowIndices, colIndices, L/(2*L-1)*thetaPartOfTerm*rowScaling(L+1))
        end
        
    end
end

% -----------------------------------------
% Add d/dzeta terms:
% -----------------------------------------

for itheta=1:Ntheta
    zetaPartOfTerm = diag(B(itheta,:))*ddzeta;
    for L=0:(Nxi-1)
        rowIndices = getIndex(itheta, 1:Nzeta, L+1, resolutionParameters);
        
        % Super-diagonal term
        if (L<Nxi-1)
            ell = L + 1;
            colIndices = getIndex(itheta,1:Nzeta,ell+1,resolutionParameters);
            addSparseBlock(rowIndices, colIndices, (L+1)/(2*L+3)*zetaPartOfTerm*rowScaling(L+1))
        end
        
        % Sub-diagonal term
        if (L>0)
            ell = L - 1;
            colIndices = getIndex(itheta,1:Nzeta,ell+1,resolutionParameters);
            addSparseBlock(rowIndices, colIndices, L/(2*L-1)*zetaPartOfTerm*rowScaling(L+1))
        end
        
    end
end

% -----------------------------------------
% Add d/dxi terms:
% -----------------------------------------

for itheta=1:Ntheta
    spatialPartOfTerm = -(iota*dBdtheta(itheta,:)+dBdzeta(itheta,:))/2;
    for L=0:(Nxi-1)
        rowIndices = getIndex(itheta, 1:Nzeta, L+1, resolutionParameters);
        
        % Super-diagonal term
        if (L<Nxi-1)
            ell = L + 1;
            colIndices = getIndex(itheta,1:Nzeta,ell+1,resolutionParameters);
            addToSparse(rowIndices, colIndices, (L+1)*(L+2)/(2*L+3)*spatialPartOfTerm*rowScaling(L+1))
        end
        
        % Sub-diagonal term
        if (L>0)
            ell = L - 1;
            colIndices = getIndex(itheta,1:Nzeta,ell+1,resolutionParameters);
            addToSparse(rowIndices, colIndices, (-L)*(L-1)/(2*L-1)*spatialPartOfTerm*rowScaling(L+1))
        end
        
    end
end

% -----------------------------------------
% Add the diffusion (collision) term:
% -----------------------------------------

L = (1:Nxi)-1;
for itheta=1:Ntheta
    for izeta=1:Nzeta
        indices = getIndex(itheta,izeta,L+1,resolutionParameters);
        addToSparse(indices, indices, nu/2*L.*(L+1).*rowScaling)
    end
end


% -----------------------------------------
% Add the extra constraint:
% -----------------------------------------

if resolutionParameters.includeConstraint
    rowIndex = matrixSize;
    L = 0;
    for itheta=1:Ntheta
        colIndices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
        addSparseBlock(rowIndex, colIndices, (1e0)*thetaWeights(itheta) * (zetaWeights') ./ (B(itheta,:) .^ 2))
    end
end

% -----------------------------------------
% Add the "source" lambda:
% -----------------------------------------

if resolutionParameters.includeConstraint
    colIndex = matrixSize;
    L = 0;
    for itheta=1:Ntheta
        rowIndices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
        addSparseBlock(rowIndices, colIndex, (1e0)*ones(Nzeta,1)*rowScaling(1))
    end
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

rhs = zeros(matrixSize,1);
spatialPart = (1./B) .* (geometryParameters.G * dBdtheta - geometryParameters.I * dBdzeta);
for itheta=1:Ntheta
    L=0;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    rhs(indices) = spatialPart(itheta,:) * (4/3) * rowScaling(L+1);
    
    L=2;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    rhs(indices) = spatialPart(itheta,:) * (2/3) * rowScaling(L+1);
end


% -----------------------------------------
% Done. Return a structure with everything important:
% -----------------------------------------

returnStruct = struct(...
    'matrixSize',matrixSize,...
    'theta',theta,...
    'zeta',zeta,...
    'theta2D',theta2D,...
    'zeta2D',zeta2D,...
    'thetaWeights',thetaWeights,...
    'zetaWeights',zetaWeights,...
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
