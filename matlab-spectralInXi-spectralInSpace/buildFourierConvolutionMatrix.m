function matrix = buildFourierConvolutionMatrix(ms, ns, B)

NN = numel(ms);

matrixSize = 2*NN-1;
assert(numel(B)==matrixSize);
matrix = zeros(matrixSize);

B = B*0.5;

for imn = 1:numel(ms)
    mi = ms(imn);
    ni = ns(imn);
    for jmn = 1:numel(ms)
        mj = ms(jmn);
        nj = ns(jmn);
        
        % If available, add the mode at m = mi + mj, n = ni + nj
        m = mi+mj;
        n = ni+nj;
        sign = 1;
        if m==0 && n<0
            n=-n;
            sign=-1;
        end
        mask = (ms==m & ns==n);
        index = find(mask);
        if numel(index)>1
            error('Should not get here!')
        elseif numel(index)==1
            % The sum mode is present.
            % i = B, j = f
            % Cos terms:
            matrix(index,   jmn)      = matrix(index,   jmn)      + B(imn);
            matrix(index,   jmn+NN-1) = matrix(index,   jmn+NN-1) - B(imn+NN-1);
            if m~=0 || n~=0  % There is no sin(0) term.
                matrix(index+NN-1,jmn)      = matrix(index+NN-1,jmn)      + sign*B(imn+NN-1);
                matrix(index+NN-1,jmn+NN-1) = matrix(index+NN-1,jmn+NN-1) + sign*B(imn);
            end
        else
            % The sum mode is outside the range we keep.
        end
        
        % If available, add the mode at m = mi - mj, n = ni - nj
        m = mi-mj;
        n = ni-nj;
        sign = 1;
        if m<0 || (m==0 && n<0)
            m=-m;
            n=-n;
            sign=-1;
        end
        mask = (ms==m & ns==n);
        index = find(mask);
        if numel(index)>1
            error('Should not get here!')
        elseif numel(index)==1
            % The difference mode is present.
            % i = B, j = f
            % Cos terms:
            matrix(index,   jmn)      = matrix(index,   jmn)      + B(imn);
            matrix(index,   jmn+NN-1) = matrix(index,   jmn+NN-1) + B(imn+NN-1);
            if m~=0 || n~=0  % There is no sin(0) term.
                matrix(index+NN-1,jmn)      = matrix(index+NN-1,jmn)      + sign*B(imn+NN-1);
                matrix(index+NN-1,jmn+NN-1) = matrix(index+NN-1,jmn+NN-1) - sign*B(imn);
            end
        else
            % The difference mode is outside the range we keep.
        end
        
        
    end
end


end