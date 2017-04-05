function matrix = mmc_populate_matrix(which_matrix, level)
global levels Nthetas Nzetas Nxis constraint_option
global G I iota nu E FSAB2 VPrime

fprintf('Entering populate_matrix: which_matrix=%d, level=%d\n',which_matrix, level)

start_time = tic;

Ntheta = Nthetas(level);
Nzeta = Nzetas(level);
Nxi = Nxis(level);
matrixSize = levels(level).matrixSize;
resolutionParameters = struct('Ntheta',Ntheta,'Nzeta',Nzeta,'Nxi',Nxi,'matrixSize',matrixSize);

switch which_matrix
    case 0
        % low order matrix
        ddtheta_plus = levels(level).ddtheta_plus_preconditioner;
        ddtheta_minus = levels(level).ddtheta_minus_preconditioner;
        ddzeta_plus = levels(level).ddzeta_plus_preconditioner;
        ddzeta_minus = levels(level).ddzeta_minus_preconditioner;
        ddxi_plus = levels(level).ddxi_plus_preconditioner;
        ddxi_minus = levels(level).ddxi_minus_preconditioner;
        pitch_angle_scattering_operator = levels(level).pitch_angle_scattering_operator_preconditioner;
    case 1
        % high order matrix
        ddtheta_plus = levels(level).ddtheta_plus;
        ddtheta_minus = levels(level).ddtheta_minus;
        ddzeta_plus = levels(level).ddzeta_plus;
        ddzeta_minus = levels(level).ddzeta_minus;
        ddxi_plus = levels(level).ddxi_plus;
        ddxi_minus = levels(level).ddxi_minus;
        pitch_angle_scattering_operator = levels(level).pitch_angle_scattering_operator;
    otherwise
        error('Invalid which_matrix')
end
B = levels(level).B;
dBdtheta = levels(level).dBdtheta;
dBdzeta = levels(level).dBdzeta;
xi = levels(level).xi;

estimated_nnz = ...
    nnz(ddtheta_plus)*Nzeta*Nxi ... % ddtheta term
    + nnz(ddzeta_plus)*Ntheta*Nxi ...  % ddzeta term
    + Ntheta*Nzeta*nnz(ddxi_plus) ... % ddxi term
    + Ntheta*Nzeta*nnz(pitch_angle_scattering_operator) ...   % collision term
    + Ntheta*Nzeta*Nxi ... % constraint term
    + Ntheta*Nzeta*Nxi; % source term

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
% Add d/dzeta terms:
% -----------------------------------------

for ixi = 1:Nxi
    for itheta=1:Ntheta
        for izetaRow=1:Nzeta
            BB = B(itheta,izetaRow);
            factor = xi(ixi)*BB - iota*E*BB*BB*I/(G*FSAB2);
            if factor > 0
                stuffToAdd = factor * ddzeta_plus(izetaRow,:);
            else
                stuffToAdd = factor * ddzeta_minus(izetaRow,:);
            end
            rowIndex = mmc_get_index(itheta, izetaRow, ixi, resolutionParameters);
            colIndices = mmc_get_index(itheta, 1:Nzeta, ixi, resolutionParameters);
            addSparseBlock(rowIndex, colIndices, stuffToAdd)
        end
    end
end

% -----------------------------------------
% Add d/dtheta terms:
% -----------------------------------------


for ixi = 1:Nxi
    for ithetaRow=1:Ntheta
        for izeta=1:Nzeta
            BB = B(ithetaRow,izeta);
            factor = iota*xi(ixi)*BB + iota*E*BB*BB/(FSAB2);
            if factor > 0
                stuffToAdd = factor * ddtheta_plus(ithetaRow,:);
            else
                stuffToAdd = factor * ddtheta_minus(ithetaRow,:);
            end
            rowIndex = mmc_get_index(ithetaRow, izeta, ixi, resolutionParameters);
            colIndices = mmc_get_index(1:Ntheta, izeta, ixi, resolutionParameters);
            addSparseBlock(rowIndex, colIndices, stuffToAdd)
        end
    end
end

% -----------------------------------------
% Add mirror term:
% -----------------------------------------

xiPart_plus = diag(1-xi.^2)*ddxi_plus;
xiPart_minus = diag(1-xi.^2)*ddxi_minus;
spatialPart = -(iota*dBdtheta+dBdzeta)/2;
for itheta=1:Ntheta
    for izeta = 1:Nzeta
        if spatialPart(itheta,izeta)>0
            xiPart_to_use = xiPart_plus;
        else
            xiPart_to_use = xiPart_minus;
        end
        indices = mmc_get_index(itheta, izeta, 1:Nxi, resolutionParameters);
        addSparseBlock(indices, indices, spatialPart(itheta,izeta)*xiPart_to_use)
    end
end

% -----------------------------------------
% Add the diffusion (collision) term:
% -----------------------------------------

for itheta=1:Ntheta
    for izeta=1:Nzeta
        indices = mmc_get_index(itheta,izeta,1:Nxi,resolutionParameters);
        addSparseBlock(indices, indices, -nu*pitch_angle_scattering_operator)
    end
end


% -----------------------------------------
% Add the extra constraint:
% -----------------------------------------

if constraint_option==1
    stuffToAdd = levels(level).xiWeights(:)';
    rowIndex = matrixSize;
    for itheta=1:Ntheta
        for izeta=1:Nzeta
            colIndices = mmc_get_index(itheta,izeta,1:Nxi,resolutionParameters);
            addSparseBlock(rowIndex, colIndices, stuffToAdd*(levels(level).thetaWeights(itheta) * levels(level).zetaWeights(izeta) / (VPrime * B(itheta,izeta)^ 2)))
        end
    end
end

% -----------------------------------------
% Add the "source" lambda:
% -----------------------------------------

if constraint_option == 1
    colIndex = matrixSize;
    for itheta=1:Ntheta
        for izeta = 1:Nzeta
            rowIndices = mmc_get_index(itheta,izeta,1:Nxi,resolutionParameters);
            addSparseBlock(rowIndices, colIndex, ones(Nxi,1))
        end
    end
end

% -----------------------------------------
% Finalize the matrix
% -----------------------------------------
matrix = createSparse();

fprintf('matrixSize: %g,  nnz: %g,  estimated nnz: %g,  sparsity: %g\n',...
    matrixSize,nnz(matrix),estimated_nnz, nnz(matrix)/(matrixSize^2))

fprintf('Took %g sec.\n',toc(start_time))
    
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
