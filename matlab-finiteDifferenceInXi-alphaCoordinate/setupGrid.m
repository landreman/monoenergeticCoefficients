function [x, D, weights] = setupGrid(N,xMin,xMax)

% This function sets up the grid and finite difference differentiation
% matrix for a periodic coordinate.

if N<5
    error('N must be at least 5.')
end

x = linspace(xMin, xMax, N+1)';
x = x(:); % Turn it into a column vector.
x(end)=[]; % Discard last point.
dx = x(2)-x(1);

weights = ones(size(x))*dx;

% Build the finite-difference differentiation matrix:
D=(-diag((1/6)*ones(N-2,1),2) ...
    + diag((4/3)*ones(N-1,1),1) ...
    - diag((4/3)*ones(N-1,1),-1) ...
    + diag((1/6)*ones(N-2,1),-2))/(2*dx);

% Due to periodicity, the corners need to be populated as well:
D(1, end) = -(4/3)/(2*dx);
D(1, end-1) = (1/6)/(2*dx);
D(2, end) = (1/6)/(2*dx);

D(end, 1) = (4/3)/(2*dx);
D(end, 2) = -(1/6)/(2*dx);
D(end-1, 1) = -(1/6)/(2*dx);

end
