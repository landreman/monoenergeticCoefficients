function returnStruct = assembleMatrix(resolutionParameters, nu, E, geometryParameters,discretizationParameters)

Nalpha = resolutionParameters.Nalpha;
Nzeta = resolutionParameters.Nzeta;
Nxi = resolutionParameters.Nxi;
buffer_zeta_points_on_each_side = discretizationParameters.buffer_zeta_points_on_each_side;

matrixSize = Nalpha*Nzeta*Nxi;
if resolutionParameters.includeConstraint
    matrixSize = matrixSize + 1;
end

% Build alpha grid:
%[alpha, ddalpha, alphaWeights] = setupGrid(Nalpha,0,2*pi);
call_uniform_diff_matrices = true;
switch abs(discretizationParameters.alpha_derivative_option)
    case 0
        fprintf('df/dalpha derivative is dropped.\n')
        call_uniform_diff_matrices = false;
        [alpha, ~, alphaWeights] = setupGrid(Nalpha,0,2*pi);
        ddalpha_plus = zeros(Nalpha);
        ddalpha_minus = zeros(Nalpha);
    case 2
        fprintf('df/dalpha derivative: centered differences, 1 point on each side.\n')
        derivative_option_plus = 0;
        derivative_option_minus = derivative_option_plus;
    case 3
        fprintf('df/dalpha derivative: centered differences, 2 points on each side.\n')
        derivative_option_plus = 10;
        derivative_option_minus = derivative_option_plus;
    case 4
        fprintf('df/dalpha derivative: upwinded differences, 0 points on one side, 1 point on the other side.\n')
        derivative_option_plus  = 30;
        derivative_option_minus = 40;
    case 5
        fprintf('df/dalpha derivative: upwinded differences, 0 points on one side, 2 points on the other side.\n')
        derivative_option_plus  = 50;
        derivative_option_minus = 60;
    case 6
        fprintf('df/dalpha derivative: upwinded differences, 1 points on one side, 2 points on the other side.\n')
        derivative_option_plus  = 80;
        derivative_option_minus = 90;
    case 7
        fprintf('df/dalpha derivative: upwinded differences, 1 point on one side, 3 points on the other side.\n')
        derivative_option_plus  = 100;
        derivative_option_minus = 110;
    case 8
        fprintf('df/dalpha derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
        derivative_option_plus  = 120;
        derivative_option_minus = 130;
end
if call_uniform_diff_matrices
    quadrature_option = 0;
    [alpha, alphaWeights, ddalpha_plus, ~]  = uniformDiffMatrices(Nalpha, ...
        0, 2*pi, derivative_option_plus, quadrature_option);
    [~, ~, ddalpha_minus, ~] = uniformDiffMatrices(Nalpha, ...
        0, 2*pi, derivative_option_minus, quadrature_option);
end
if discretizationParameters.alpha_derivative_option<0
    ddalpha_plus = diag(diag(ddalpha_plus));
    ddalpha_minus = diag(diag(ddalpha_minus));
    fprintf('  But only the diagonal is kept.\n')
end



% Build zeta grid:
zetaMax = 2*pi/geometryParameters.Nperiods;

switch abs(discretizationParameters.zeta_derivative_option)
    case 2
        fprintf('df/dzeta derivative: centered differences, 1 point on each side.\n')
        derivative_option_plus = 2;
        derivative_option_minus = derivative_option_plus;
    case 3
        fprintf('df/dzeta derivative: centered differences, 2 points on each side.\n')
        derivative_option_plus = 12;
        derivative_option_minus = derivative_option_plus;
    case 4
        fprintf('df/dzeta derivative: upwinded differences, 0 points on one side, 1 point on the other side.\n')
        derivative_option_plus  = 32;
        derivative_option_minus = 42;
    case 5
        fprintf('df/dzeta derivative: upwinded differences, 0 points on one side, 2 points on the other side.\n')
        derivative_option_plus  = 52;
        derivative_option_minus = 62;
    case 6
        fprintf('df/dzeta derivative: upwinded differences, 1 points on one side, 2 points on the other side.\n')
        derivative_option_plus  = 82;
        derivative_option_minus = 92;
    case 7
        fprintf('df/dzeta derivative: upwinded differences, 1 point on one side, 3 points on the other side.\n')
        derivative_option_plus  = 102;
        derivative_option_minus = 112;
    case 8
        fprintf('df/dzeta derivative: upwinded differences, 2 points on one side, 3 points on the other side.\n')
        derivative_option_plus  = 122;
        derivative_option_minus = 132;
end
fprintf('buffer_zeta_points_on_each_side: %d\n',buffer_zeta_points_on_each_side)
Delta = (2*pi)/(geometryParameters.Nperiods*(Nzeta-2*buffer_zeta_points_on_each_side));
quadrature_option = 0;
[zeta, ~, ddzeta_plus, ~]  = uniformDiffMatrices(Nzeta, ...
    -buffer_zeta_points_on_each_side*Delta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta, derivative_option_plus, quadrature_option);
[~   , ~, ddzeta_minus, ~] = uniformDiffMatrices(Nzeta, ...
    -buffer_zeta_points_on_each_side*Delta, zetaMax+(buffer_zeta_points_on_each_side-1)*Delta, derivative_option_minus, quadrature_option);
assert(abs(zeta(2)-zeta(1)-Delta)<1e-12)

if discretizationParameters.zeta_derivative_option<0
    fprintf('  But only the diagonal is kept.\n')
    ddzeta_plus = diag(diag(ddzeta_plus));
    ddzeta_minus = diag(diag(ddzeta_minus));
end

ddzeta_sum        = 0.5*(ddzeta_plus+ddzeta_minus);
ddzeta_difference = 0.5*(ddzeta_plus-ddzeta_minus);

zetaWeights=ones(size(zeta));
zetaWeights(1:buffer_zeta_points_on_each_side)         = 0;
zetaWeights(end-buffer_zeta_points_on_each_side+1:end) = 0;
zetaWeights = zetaWeights * Delta * geometryParameters.Nperiods;
assert(abs(sum(zetaWeights)-2*pi) < 1e-12)

zeta_to_impose_DKE = (buffer_zeta_points_on_each_side+1):(Nzeta-buffer_zeta_points_on_each_side);

[zeta2D, alpha2D] = meshgrid(zeta,alpha);
theta2D = alpha2D + geometryParameters.iota*zeta2D;

% Initialize the arrays for the magnetic field strength B:
iota = geometryParameters.iota;
[B, dBdtheta, dBdzeta] = geometry(theta2D, zeta2D, geometryParameters);

%{
estimated_nnz = ...
    Nalpha*Nzeta*Nxi*4*2 ... % ddalpha term (4 off-diagonals in alpha, 2 off-diagonals in xi)
    + Nalpha*Nzeta*Nxi*4*2 ...  % ddzeta term (4 off-diagonals in zeta, 2 off-diagonals in xi)
    + Nalpha*Nzeta*Nxi*2 ... % ddxi term (2 off-diagonals in xi)
    + Nalpha*Nzeta*(Nxi-1) ...   % collision term (diagonal, 0 when L=0)
    + Nalpha*Nzeta ... % constraint term
    + Nalpha*Nzeta; % source term
%}
estimated_nnz = ...
    (1+discretizationParameters.alpha_interpolation_stencil)*Nalpha*(2*buffer_zeta_points_on_each_side)*Nxi ... % alpha interpolation
    + 2*Nalpha*Nxi*nnz(ddzeta_sum) ...  % ddzeta_sum term (2 off-diagonals in xi)
    + 3*Nalpha*Nxi*nnz(ddzeta_difference) ...  % ddzeta_difference term (3 off-diagonals in xi)
    + Nalpha*Nzeta*Nxi*2 ... % ddxi term (2 off-diagonals in xi)
    + Nalpha*Nzeta*(Nxi-1) ...   % collision term (diagonal, 0 when L=0)
    + Nalpha*Nzeta ... % constraint term
    + Nalpha*Nzeta; % source term

% Initialize some arrays that are used for building the sparse matrix:
sparseCreatorIndex=1;
sparseCreator_i=0;
sparseCreator_j=0;
sparseCreator_s=0;
resetSparseCreator()

VPrime = alphaWeights' * (1./B.^2) * zetaWeights;
FSAB2 = (alphaWeights' * ones(size(B)) * zetaWeights) / VPrime;

if discretizationParameters.piecewise_upwinding
    fprintf('Computing Legendre integrals.\n')
    LegendreIntegrals0 = zeros(Nxi,1);
    LegendreIntegrals1 = zeros(Nxi,1);
    for L=0:(Nxi-1)
        LegendreIntegrals0(L+1) = integral(@integrand0,0,1);
        LegendreIntegrals1(L+1) = integral(@integrand1,0,1);
    end
    
    Ls=(1:Nxi)-1;
    LegendreIntegrals2 = (Ls+1)./((2*Ls+1).*(2*Ls+3));
    LegendreIntegrals1-LegendreIntegrals2(:)
end

    function zz=integrand0(xi)
        P = legendre(L,xi);
        P = P(1,:);
        P = reshape(P,size(xi));
        zz = xi.*P.*P;
    end

    function zz=integrand1(xi)
        P = legendre(L,xi);
        P = P(1,:);
        Q = legendre(L+1,xi);
        Q = Q(1,:);
        PQ = reshape(P.*Q,size(xi));
        zz = xi.*PQ;
    end

% ***************************************************************************
% ***************************************************************************
% Now populate the matrix.
% ***************************************************************************
% ***************************************************************************

Ls = (1:Nxi)-1;
%rowScaling = (2*Ls+1) .^ (1.5);
rowScaling = ones(size(Ls));


% -----------------------------------------
% Add d/dzeta terms:
% -----------------------------------------
for ialpha=1:Nalpha
    zetaPartOfTerm_full = diag(B(ialpha,:))*ddzeta_sum;
    zetaPartOfTerm_sum = zetaPartOfTerm_full(zeta_to_impose_DKE,:);
    
    zetaPartOfTerm_full = diag(B(ialpha,:))*ddzeta_difference;
    zetaPartOfTerm_difference_diagInXi = zetaPartOfTerm_full(zeta_to_impose_DKE,:);
    if discretizationParameters.xi_pentadiagonal_terms_option==2
        zetaPartOfTerm_full = diag(diag(zetaPartOfTerm_full));
    end
    zetaPartOfTerm_difference_pentadiagInXi = zetaPartOfTerm_full(zeta_to_impose_DKE,:);
    
    for L=0:(Nxi-1)
        rowIndices = getIndex(ialpha, zeta_to_impose_DKE, L+1, resolutionParameters);
        
        if ~ discretizationParameters.piecewise_upwinding
            % Diagonal term
            ell = L;
            colIndices = getIndex(ialpha,1:Nzeta,ell+1,resolutionParameters);
            addSparseBlock(rowIndices, colIndices, (2*L*L+2*L-1)/((2*L+3)*(2*L-1))*zetaPartOfTerm_difference_diagInXi*rowScaling(L+1))
            
            if discretizationParameters.include_streaming_xi_tridiagonal_terms
                % Super-diagonal term
                if (L<Nxi-1)
                    ell = L + 1;
                    colIndices = getIndex(ialpha,1:Nzeta,ell+1,resolutionParameters);
                    addSparseBlock(rowIndices, colIndices, (L+1)/(2*L+3)*zetaPartOfTerm_sum*rowScaling(L+1))
                end
                
                % Sub-diagonal term
                if (L>0)
                    ell = L - 1;
                    colIndices = getIndex(ialpha,1:Nzeta,ell+1,resolutionParameters);
                    addSparseBlock(rowIndices, colIndices, L/(2*L-1)*zetaPartOfTerm_sum*rowScaling(L+1))
                end
            end
            
            if discretizationParameters.xi_pentadiagonal_terms_option ~= 0
                % Super-super-diagonal term
                if (L<Nxi-2)
                    ell = L + 2;
                    colIndices = getIndex(ialpha,1:Nzeta,ell+1,resolutionParameters);
                    addSparseBlock(rowIndices, colIndices, (L+2)*(L+1)/((2*L+5)*(2*L+3))*zetaPartOfTerm_difference_pentadiagInXi*rowScaling(L+1))
                end
                
                % Sub-sub-diagonal term
                if (L>1)
                    ell = L - 2;
                    colIndices = getIndex(ialpha,1:Nzeta,ell+1,resolutionParameters);
                    addSparseBlock(rowIndices, colIndices, (L-1)*L/((2*L-3)*(2*L-1))*zetaPartOfTerm_difference_pentadiagInXi*rowScaling(L+1))
                end
            end
        else
            % DO use the piecewise-interval upwinding
            
            % Diagonal term
            ell = L;
            colIndices = getIndex(ialpha,1:Nzeta,ell+1,resolutionParameters);
            addSparseBlock(rowIndices, colIndices, LegendreIntegrals0(L+1)*(2*L+1)/(2)*zetaPartOfTerm_difference_diagInXi*rowScaling(L+1))
        end
    end
end

% -----------------------------------------
% Add d/alpha term:
% -----------------------------------------
factor = (geometryParameters.iota*E/FSAB2)*(1 + geometryParameters.iota*geometryParameters.I/geometryParameters.G);
if factor>0
    ddalpha_to_use = ddalpha_plus;
else
    ddalpha_to_use = ddalpha_minus;
end
for izeta = zeta_to_impose_DKE
    stuffToAdd = factor*diag(B(:,izeta).^2)*ddalpha_to_use;
    for L=0:(Nxi-1)
        indices = getIndex(1:Nalpha, izeta, L+1, resolutionParameters);
        addSparseBlock(indices, indices, stuffToAdd*rowScaling(L+1))
    end
end

% -----------------------------------------
% Add periodicity constraints:
% -----------------------------------------

% First handle the points needed by the DKE to the left of zeta=0:
theta_left = alpha;
theta_right = alpha - iota * zetaMax;
switch discretizationParameters.alpha_interpolation_stencil
    case 100
        interpolationMatrix = eye(Nalpha);
        %interpolationMatrix = zeros(Nalpha);
    case 0
        interpolationMatrix = m20130226_06_FourierSpectralInterpolationMatrix(Nalpha, theta_right);
    case {1,2,3,4,5}
        interpolationMatrix = m20160925_01_periodicInterpolation(theta_left, theta_right, 2*pi, discretizationParameters.alpha_interpolation_stencil);
    otherwise
        error('invalid alpha_interpolation_stencil')
end
izetas = 1:buffer_zeta_points_on_each_side;
izeta_shift = Nzeta-2*buffer_zeta_points_on_each_side;
if discretizationParameters.alpha_interpolation_stencil<3
    assignin('base','alpha_interpolation_matrix_left',interpolationMatrix)
end
for ixi = 1:Nxi
    if ixi>discretizationParameters.preconditioner_alpha_min_L
    %if false
        %interpolationMatrixToUse = zeros(Nalpha);
        interpolationMatrixToUse = eye(Nalpha);
    else
        interpolationMatrixToUse = interpolationMatrix;
    end
    for izeta = izetas
        % Add 1's along the diagonal
        rowIndices = getIndex(1:Nalpha, izeta, ixi, resolutionParameters);
        addToSparse(rowIndices, rowIndices, ones(size(rowIndices)));
        %if ~ preconditioner
            % Add interpolation matrix
            colIndices = getIndex(1:Nalpha, izeta+izeta_shift, ixi, resolutionParameters);
            addSparseBlock(rowIndices,colIndices,-interpolationMatrixToUse);
        %end
    end
end

% Now handle the points needed by the DKE to the right of zeta=zetaMax-Delta:
theta_left = alpha;
theta_right = alpha + iota * zetaMax;
switch discretizationParameters.alpha_interpolation_stencil
    case 100
        interpolationMatrix = eye(Nalpha);
        %interpolationMatrix = zeros(Nalpha);
        %interpolationMatrix = m20160925_01_periodicInterpolation(theta_left, theta_right, 2*pi, 2);
    case 0
        interpolationMatrix = m20130226_06_FourierSpectralInterpolationMatrix(Nalpha, theta_right);
    case {1,2,3,4,5}
        interpolationMatrix = m20160925_01_periodicInterpolation(theta_left, theta_right, 2*pi, discretizationParameters.alpha_interpolation_stencil);
    otherwise
        error('invalid alpha_interpolation_stencil')
end
if discretizationParameters.alpha_interpolation_stencil<3
    assignin('base','alpha_interpolation_matrix_right',interpolationMatrix)
end
izetas = (Nzeta-buffer_zeta_points_on_each_side+1):Nzeta;

for ixi = 1:Nxi
    if ixi>discretizationParameters.preconditioner_alpha_min_L
    %if ixi>preconditioner_min_L && preconditioner
    %if false
        %interpolationMatrixToUse = zeros(Nalpha);
        interpolationMatrixToUse = eye(Nalpha);
    else
        interpolationMatrixToUse = interpolationMatrix;
    end
    for izeta = izetas
        % Add 1's along the diagonal
        rowIndices = getIndex(1:Nalpha, izeta, ixi, resolutionParameters);
        addToSparse(rowIndices, rowIndices, ones(size(rowIndices)));
        %if ~ preconditioner
            % Add interpolation matrix
            colIndices = getIndex(1:Nalpha, izeta-izeta_shift, ixi, resolutionParameters);
            addSparseBlock(rowIndices,colIndices,-interpolationMatrixToUse);
        %end
    end
end

% -----------------------------------------
% Add d/dxi terms:
% -----------------------------------------

for ialpha=1:Nalpha
    spatialPartOfTerm = -(iota*dBdtheta(ialpha,zeta_to_impose_DKE)+dBdzeta(ialpha,zeta_to_impose_DKE))/2;
    for L=0:(Nxi-1)
        rowIndices = getIndex(ialpha, zeta_to_impose_DKE, L+1, resolutionParameters);
        
        % Super-diagonal term
        if (L<Nxi-1)
            ell = L + 1;
            colIndices = getIndex(ialpha,zeta_to_impose_DKE,ell+1,resolutionParameters);
            addToSparse(rowIndices, colIndices, (L+1)*(L+2)/(2*L+3)*spatialPartOfTerm*rowScaling(L+1))
        end
        
        % Sub-diagonal term
        if (L>0)
            ell = L - 1;
            colIndices = getIndex(ialpha,zeta_to_impose_DKE,ell+1,resolutionParameters);
            addToSparse(rowIndices, colIndices, (-L)*(L-1)/(2*L-1)*spatialPartOfTerm*rowScaling(L+1))
        end
        
    end
end

% -----------------------------------------
% Add the diffusion (collision) term:
% -----------------------------------------

L = (1:Nxi)-1;
stuffToAdd = nu/2*L.*(L+1).*rowScaling;
%if discretizationParameters.alpha_interpolation_stencil==100 || discretizationParameters.alpha_interpolation_stencil==1
if discretizationParameters.alpha_interpolation_stencil==100
    stuffToAdd(1)=stuffToAdd(2);
end
for ialpha=1:Nalpha
    for izeta=zeta_to_impose_DKE
        indices = getIndex(ialpha,izeta,L+1,resolutionParameters);
        %addToSparse(indices, indices, nu/2*L.*(L+1).*rowScaling)
        addToSparse(indices, indices, stuffToAdd)
    end
end


% -----------------------------------------
% Add the extra constraint:
% -----------------------------------------

if resolutionParameters.includeConstraint
    rowIndex = matrixSize;
    L = 0;
    for ialpha=1:Nalpha
        colIndices = getIndex(ialpha,1:Nzeta,L+1,resolutionParameters);
        addSparseBlock(rowIndex, colIndices, (1e0)*alphaWeights(ialpha) * (zetaWeights') ./ (B(ialpha,:) .^ 2))
    end
end

% -----------------------------------------
% Add the "source" lambda:
% -----------------------------------------

if resolutionParameters.includeConstraint
    colIndex = matrixSize;
    L = 0;
    for ialpha=1:Nalpha
        rowIndices = getIndex(ialpha,zeta_to_impose_DKE,L+1,resolutionParameters);
        addSparseBlock(rowIndices, colIndex, (1e0)*ones(numel(zeta_to_impose_DKE),1)*rowScaling(1))
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
for ialpha=1:Nalpha
    L=0;
    indices = getIndex(ialpha,zeta_to_impose_DKE,L+1,resolutionParameters);
    rhs(indices) = spatialPart(ialpha,zeta_to_impose_DKE) * (4/3) * rowScaling(L+1);
    
    L=2;
    indices = getIndex(ialpha,zeta_to_impose_DKE,L+1,resolutionParameters);
    rhs(indices) = spatialPart(ialpha,zeta_to_impose_DKE) * (2/3) * rowScaling(L+1);
end


% -----------------------------------------
% Done. Return a structure with everything important:
% -----------------------------------------

returnStruct = struct(...
    'matrixSize',matrixSize,...
    'alpha',alpha,...
    'zeta',zeta,...
    'alpha2D',alpha2D,...
    'theta2D',theta2D,...
    'zeta2D',zeta2D,...
    'alphaWeights',alphaWeights,...
    'zetaWeights',zetaWeights,...
    'B',B,...
    'dBdtheta',dBdtheta,...
    'dBdzeta',dBdzeta,...
    'ddzeta_plus',ddzeta_plus,...
    'ddzeta_minus',ddzeta_minus,...
    'ddzeta_sum',ddzeta_sum,...
    'ddzeta_difference',ddzeta_difference,...
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
