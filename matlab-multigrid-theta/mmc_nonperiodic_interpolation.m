function m = nonperiodicInterpolation(x,y)
% x = old grid
% y = new grid

% Verify the old grid x is in increasing order, with no repeated entried:
differences = x(2:end) - x(1:end-1);
assert(all(differences>0))

m = zeros(numel(y), numel(x));

for i=1:numel(y)
    yy = y(i);
    if yy <= x(1)
        m(i,1) = 1;
    elseif yy >= x(end)
        m(i,end) = 1;
    else
        index = find(x>yy,1);
        if numel(index)>1
            index
            error('numel(index)>1')
        elseif numel(index)==0
            error('numel(index) == 0')
        elseif index==1
            % The 'if yy <= x(1)' case above should prevent this from
            % happening.
            error('Should not get here!')
        end
        indicesToUse = [-1,0] + index;
        
        x0 = x(indicesToUse(1));
        x1 = x(indicesToUse(2));
        assert(yy >= x0)
        assert(yy <= x1)
        
        m(i,indicesToUse(1)) = (yy-x1)/( (x0-x1) );
        m(i,indicesToUse(2)) = (yy-x0)/( (x1-x0) );

    end
end

end