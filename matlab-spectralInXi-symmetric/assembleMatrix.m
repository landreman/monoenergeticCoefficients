function returnStruct = assembleMatrix(resolutionParameters, nu_prime, E, geometryParameters, derivative_option)

Ntheta = resolutionParameters.Ntheta;
Nzeta = resolutionParameters.Nzeta;
Nxi = resolutionParameters.Nxi;
includeConstraints = resolutionParameters.includeConstraints;

matrixSize = Ntheta*Nzeta*Nxi + Ntheta*Nzeta;
if includeConstraints
    matrixSize = matrixSize + 2;
end

nu_hat = nu_prime / (geometryParameters.G + geometryParameters.iota * geometryParameters.I);

% Build theta grid:
[theta, ddtheta, thetaWeights] = setupGrid(Ntheta,0,2*pi,derivative_option);

% Build zeta grid:
if Nzeta==1
    zeta=0;
    ddzeta=0;
    zetaWeights=2*pi;
else
    zetaMax = 2*pi/geometryParameters.Nperiods;
    [zeta, ddzeta, zetaWeights] = setupGrid(Nzeta,0,zetaMax,derivative_option);
    zetaWeights = zetaWeights * geometryParameters.Nperiods;
end

[zeta2D, theta2D] = meshgrid(zeta,theta);

% Integration weight should be independent of theta and zetao:
thetaWeight = thetaWeights(1);
zetaWeight  = zetaWeights(1);

% Initialize the arrays for the magnetic field strength B:
iota = geometryParameters.iota;
[B, dBdtheta, dBdzeta] = geometry(theta2D, zeta2D, geometryParameters);

sqrt_g = (geometryParameters.G + geometryParameters.iota*geometryParameters.I) ./ (B.*B);

VPrime = thetaWeights' * sqrt_g * zetaWeights;
FSAB2 = thetaWeights' * (sqrt_g.*B.*B) * zetaWeights / VPrime;
fprintf('FSAB2: %g\n',FSAB2)

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
% Set scaling
% ***************************************************************************
% ***************************************************************************

Ls = (1:Nxi)-1;
L_scaling = ones(size(Ls));
%L_scaling = sqrt((2*Ls+1)/2);

%assert(abs(L_scaling(2)^2 - 3*L_scaling(1)^2) < 1e-12)

% ***************************************************************************
% ***************************************************************************
% Create a diagonal matrix that has the effect of applying < \int dxi (...) >
% ***************************************************************************
% ***************************************************************************

L = (1:Nxi)-1;
% In the line below, L_scaling^2 converts p's to P's, and 2/(2L+1) is the result of integrating two P's.
L_factors = L_scaling.*L_scaling.*(2./(2*L+1.0));
factor = (thetaWeight*zetaWeight/VPrime)*sqrt_g; % Including the weights for < ... >
for itheta=1:Ntheta
    for izeta=1:Nzeta
        indices = getIndex(itheta,izeta,L+1,resolutionParameters);
        addToSparse(indices, indices, factor(itheta,izeta)*L_factors)
    end
end

weights_matrix = createSparse();

% ***************************************************************************
% ***************************************************************************
% Populate the collisionless part of the drift-kinetic operator (V)
% ***************************************************************************
% ***************************************************************************

% -----------------------------------------
% Add d/dtheta terms:
% -----------------------------------------

spatial_part_of_streaming_term = iota ./ (B .* sqrt_g);
spatial_part_of_ExB_term = (E * iota) ./ (sqrt_g * FSAB2);
for izeta=1:Nzeta
    theta_part_of_streaming_term = diag(spatial_part_of_streaming_term(:,izeta))*ddtheta;
    theta_part_of_ExB_term = diag(spatial_part_of_ExB_term(:,izeta))*ddtheta;
    for L=0:(Nxi-1)
        rowIndices = getIndex(1:Ntheta,izeta,L+1,resolutionParameters);
        
        % Diagonal term (ExB)
        ell = L;
        colIndices = rowIndices;
        addSparseBlock(rowIndices, colIndices, theta_part_of_ExB_term)
        
        % Super-diagonal term (streaming)
        if (L<Nxi-1)
            ell = L + 1;
            colIndices = getIndex(1:Ntheta,izeta,ell+1,resolutionParameters);
            addSparseBlock(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*ell/(2*ell+1.0)*theta_part_of_streaming_term)
            %addSparseBlock(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*(ell+1)/(2*ell+1.0)*theta_part_of_streaming_term)
        end
        
        % Sub-diagonal term
        if (L>0)
            ell = L - 1;
            colIndices = getIndex(1:Ntheta,izeta,ell+1,resolutionParameters);
            addSparseBlock(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*(ell+1)/(2*ell+1.0)*theta_part_of_streaming_term)
            %addSparseBlock(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*ell/(2*ell+1.0)*theta_part_of_streaming_term)
        end
        
    end
end

% -----------------------------------------
% Add d/dzeta terms:
% -----------------------------------------

spatial_part_of_streaming_term = 1 ./ (B .* sqrt_g);
spatial_part_of_ExB_term = (-iota * E * geometryParameters.I / geometryParameters.G) ./ (sqrt_g * FSAB2);
for itheta=1:Ntheta
    zeta_part_of_streaming_term = diag(spatial_part_of_streaming_term(itheta,:))*ddzeta;
    zeta_part_of_ExB_term = diag(spatial_part_of_ExB_term(itheta,:))*ddzeta;
    for L=0:(Nxi-1)
        rowIndices = getIndex(itheta, 1:Nzeta, L+1, resolutionParameters);
        
        % Diagonal term (ExB)
        ell = L;
        colIndices = rowIndices;
        addSparseBlock(rowIndices, colIndices, zeta_part_of_ExB_term)

        % Super-diagonal term (streaming)
        if (L<Nxi-1)
            ell = L + 1;
            colIndices = getIndex(itheta,1:Nzeta,ell+1,resolutionParameters);
            %addSparseBlock(rowIndices, colIndices, (L_scaling(L)/L_scaling(ell))*(L+1)/(2*L+1.0)*zeta_part_of_streaming_term)
            %addSparseBlock(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*(ell+1)/(2*ell+1.0)*zeta_part_of_streaming_term)
            addSparseBlock(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*ell/(2*ell+1.0)*zeta_part_of_streaming_term)
        end
        
        % Sub-diagonal term (streaming)
        if (L>0)
            ell = L - 1;
            colIndices = getIndex(itheta,1:Nzeta,ell+1,resolutionParameters);
            %addSparseBlock(rowIndices, colIndices, (L_scaling(L)/L_scaling(ell))*L/(2*L+1.0)*zeta_part_of_streaming_term)
            %addSparseBlock(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*ell/(2*ell+1.0)*zeta_part_of_streaming_term)
            addSparseBlock(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*(ell+1)/(2*ell+1.0)*zeta_part_of_streaming_term)
        end
        
    end
end

% -----------------------------------------
% Add d/dxi terms:
% -----------------------------------------

spatialPartOfTerm = -(iota*dBdtheta+dBdzeta) ./ (2*B.*B.*sqrt_g);
for itheta=1:Ntheta
    for L=0:(Nxi-1)
        rowIndices = getIndex(itheta, 1:Nzeta, L+1, resolutionParameters);
        
        % Super-diagonal term
        if (L<Nxi-1)
            ell = L + 1;
            colIndices = getIndex(itheta,1:Nzeta,ell+1,resolutionParameters);
            %addToSparse(rowIndices, colIndices, -(L_scaling(L)/L_scaling(ell))*L*(L+1)/(2*L+1)*spatialPartOfTerm(itheta,:))
            addToSparse(rowIndices, colIndices, (L_scaling(ell+1)/L_scaling(L+1))*ell*(ell+1)/(2*ell+1)*spatialPartOfTerm(itheta,:))
        end
        
        % Sub-diagonal term
        if (L>0)
            ell = L - 1;
            colIndices = getIndex(itheta,1:Nzeta,ell+1,resolutionParameters);
            %addToSparse(rowIndices, colIndices, (L_scaling(L)/L_scaling(ell))*L*(L+1)/(2*L+1)*spatialPartOfTerm(itheta,:))
            addToSparse(rowIndices, colIndices, -(L_scaling(ell+1)/L_scaling(L+1))*ell*(ell+1)/(2*ell+1)*spatialPartOfTerm(itheta,:))
        end
        
    end
end

% ***************************************************************************
% ***************************************************************************
% Now populate the matrix.
% ***************************************************************************
% ***************************************************************************

V = createSparse();

weights_matrix_V = weights_matrix*V;
weights_matrix_V_symm = weights_matrix_V + weights_matrix_V';

% ***************************************************************************
% Create a matrix for the 'real' collision term, including the weights for < \int dxi (...) >
% ***************************************************************************

L = (1:Nxi)-1;
% In the line below, -L*(L+1) is from the collision operator, L_scaling^2
% converts p's to P's, and 2/(2L+1) is the result of integrating two P's.
L_factors = -L.*(L+1).*L_scaling.*L_scaling.*(2./(2*L+1.0));
% This next line assumes thetaWeights and zetaWeights are constant
factor = ((nu_hat/2)*thetaWeight*zetaWeight/VPrime)*sqrt_g; % Including the weights for < \int dxi (...) >
for itheta=1:Ntheta
    for izeta=1:Nzeta
        indices = getIndex(itheta,izeta,L+1,resolutionParameters);
        addToSparse(indices, indices, factor(itheta,izeta)*L_factors)
    end
end

real_collision_operator = createSparse();

% ***************************************************************************
% Create a matrix for \hat{C}^{-1}, including the weights for < \int dxi (...) >
% This matrix gives 0 when acting on the P_0 (constant) Legendre polynomial.
% ***************************************************************************


L = (1:Nxi)-1;
% In the line below, -1/(L*(L+1)) is from the collision operator, L_scaling^2
% converts p's to P's, and 2/(2L+1) is the result of integrating two P's.
L_factors = -(2 ./ (nu_hat * L .* (L+1))) .* L_scaling.*L_scaling.*(2./(2*L+1.0)) * thetaWeight * zetaWeight / VPrime;
% The line above assumes thetaWeights and zetaWeights are constant.
L_factors(1) = 0;
for itheta=1:Ntheta
    for izeta=1:Nzeta
        indices = getIndex(itheta,izeta,L+1,resolutionParameters);
        addToSparse(indices, indices, L_factors * sqrt_g(itheta,izeta))
    end
end

CHat_inverse = createSparse();

% Build the term in the matrix corresponding to W f(0) = V \hat{C}^{-1} V f(0) - C f(0):
matrix = -(V')*CHat_inverse * V - real_collision_operator;

%{
% ***************************************************************************
% Create a matrix for C_0^{-1}, including the weights for < \int dxi (...) >
% ***************************************************************************


% In the line below, L_scaling^2 converts p's to P's, and 2/(2L+1) is the result of integrating two P's.
L=0;
L_factors = -(2/nu) * L_scaling(1)*L_scaling(1)*(2./(2*L + 1.0)) * thetaWeights(1) * zetaWeights(1);
% The line above assumes thetaWeights and zetaWeights are constant.
for itheta=1:Ntheta
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    addToSparse(indices, indices, L_factors * sqrt_g(itheta,:))
end

C0_inverse = createSparse();

V_C0inv_V = -(V')*C0_inverse * V;
%}

% ***************************************************************************
% Take the L=0 terms of V, representing mass conservation, and add them as
% extra rows to the main matrix.
% ***************************************************************************

indices_for_L0 = zeros(Ntheta*Nzeta,1);
L=0;
for itheta = 1:Ntheta
    indices_for_L0((1:Nzeta) + (itheta-1)*Nzeta) = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
end
indices_for_L0 = sort(indices_for_L0);

% Extract the V operator operating on p_0:
V_on_L0 = V(:,indices_for_L0);
stuff_to_add = - weights_matrix * V_on_L0;
matrix(:,indices_for_L0 + Ntheta*Nzeta*Nxi) = stuff_to_add;
matrix(indices_for_L0 + Ntheta*Nzeta*Nxi,:) = stuff_to_add';

% ***************************************************************************
% Handle sources and constraints
% ***************************************************************************

if includeConstraints
    % These next 2 lines could be swapped. Would it make a difference?
    rowIndex1 = matrixSize-1;
    rowIndex2 = matrixSize;
    L = 0;
    for itheta=1:Ntheta
        stuff_to_add = (thetaWeight * zetaWeight / VPrime) * sqrt_g(itheta,:);
        
        colIndices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
        addSparseBlock(rowIndex1, colIndices, stuff_to_add)
        addSparseBlock(colIndices, rowIndex1, stuff_to_add')
        
        colIndices = colIndices + Ntheta*Nzeta*Nxi;
        addSparseBlock(rowIndex2, colIndices, stuff_to_add)
        addSparseBlock(colIndices, rowIndex2, stuff_to_add')
    end
end

% -----------------------------------------
% Finalize the matrix
% ------------------------------------------
matrix = matrix + createSparse();

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

% First build R in the space of the normalized Legendre polynomials p_L:
% (little p's), NOT including the < \int dxi (...) >
rhs_of_original_dke = zeros(matrixSize,1);
spatialPart = (1./(B.*B.*B.*sqrt_g)) .* (geometryParameters.G * dBdtheta - geometryParameters.I * dBdzeta);
for itheta=1:Ntheta
    L=0;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    rhs_of_original_dke(indices) = spatialPart(itheta,:) * (4/3) / L_scaling(L+1);
    
    L=2;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    rhs_of_original_dke(indices) = spatialPart(itheta,:) * (2/3) / L_scaling(L+1);
end

% Also form a version that DOES include the < \int dxi (...) >
R_with_weights = weights_matrix * rhs_of_original_dke;

% Next, form the terms R + V CHat^{-1} R:
rhs = R_with_weights + (-V') * CHat_inverse * rhs_of_original_dke;
% In the above line, CHat_inverse contains the factors for the < \int dxi (...) >.
% We don't want the V to act on those factors, so we integrate by parts,
% which has the effect of replacing V -> -V'.

% Add R_0 to the mathcal{F}(1) equation:
L = 0;
for itheta = 1:Ntheta
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    rhs(indices + Ntheta*Nzeta*Nxi) = R_with_weights(indices);
end

%{
% Now add R to the rhs of the main equation, including < \int dxi (...) >
spatialPart = (1./(B.*B.*B.*sqrt_g)) .* (geometryParameters.G * dBdtheta - geometryParameters.I * dBdzeta) ...  % This line is for R itself
    .* sqrt_g * (thetaWeight * zetaWeight / VPrime);                                                            % This line adds the < ... >
for itheta=1:Ntheta
    L=0;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    rhs(indices) = rhs(indices) + spatialPart(itheta,:) * (4/3) / L_scaling(L+1) ... % This line is for R itself
        * L_scaling(L+1) * L_scaling(L+1) * (2/(2*L+1.0));                           % This line adds the \int dxi
    
    L=2;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters);
    rhs(indices) = rhs(indices) + spatialPart(itheta,:) * (2/3) / L_scaling(L+1) ... % This line is for R itself
        * L_scaling(L+1) * L_scaling(L+1) * (2/(2*L+1.0));                           % This line adds the \int dxi
    
    % Also add R_0 to the mathcal{F}(1) equation:
    L=0;
    indices = getIndex(itheta,1:Nzeta,L+1,resolutionParameters) + Ntheta*Nzeta*Nxi;
    rhs(indices) = rhs(indices) + spatialPart(itheta,:) * (4/3) / L_scaling(L+1) ... % This line is for R itself
        * L_scaling(L+1) * L_scaling(L+1) * (2/(2*L+1.0));                           % This line adds the \int dxi
    
end
%}


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
    'sqrt_g',sqrt_g,...
    'L_scaling',L_scaling,...
    'weights_matrix',weights_matrix,...
    'matrix',matrix,...
    'V',V,...
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
