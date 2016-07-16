function [ddtheta,ddzeta] = buildFourierDifferentiationMatrices(ms,ns)

assert(ms(1)==0)
assert(ns(1)==0)

NN = numel(ns);
matrixSize = 2*NN-1;

%{
rowIndices = (1:matrixSize)';
colIndices = [NN+(1:(NN-1))'; (1:NN)'];
values = [ms(2:end); -ms];
ddtheta = sparse(rowIndices, colIndices, values);

rowIndices = (1:matrixSize)';
colIndices = [NN+(1:(NN-1))'; (1:NN)'];
values = [-ns(2:end); ns];
ddzeta = sparse(rowIndices, colIndices, values);
%}

rowIndices = (1:matrixSize)';
colIndices = [NN-1+(1:NN)'; (2:NN)'];
values = [ms; -ms(2:end)];
ddtheta = sparse(rowIndices, colIndices, values);

rowIndices = (1:matrixSize)';
colIndices = [NN-1+(1:NN)'; (2:NN)'];
values = [-ns; ns(2:end)];
ddzeta = sparse(rowIndices, colIndices, values);

end