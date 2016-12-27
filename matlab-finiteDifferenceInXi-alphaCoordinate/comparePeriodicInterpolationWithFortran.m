stencil=3;
period=2*pi;

for N=4:7
    x = linspace(0,2*pi,N+1);
    x(end)=[];
    for shift = (-7.2):2.1:7.2
    %for shift = (-1):1
        y = x + shift;
        matrix = m20160925_01_periodicInterpolation(x,y,period,stencil);
        fprintf('N=%d, shift=%g\n',N,shift)
        matrix
    end
end