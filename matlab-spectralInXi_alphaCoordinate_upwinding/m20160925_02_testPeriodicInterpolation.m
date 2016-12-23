function m20160925_02_testPeriodicInterpolation()

stencil = 0;

    function ff = f(x)
        ff = cos(1+0.9*cos(x+1))+log(1.5+sin(2*x-1));
    end

Nalpha = 15;
shift = 2.2;

theta1 = linspace(0,2*pi,Nalpha+1)';
theta1(end)=[];
%theta1 = theta1+0.99*theta1(2);

theta2 = theta1+shift;
%theta2 = mod(theta2,2*pi);

theta2 = linspace(0,2*pi,300)' ;
%theta2 = 2*pi-linspace(0,0.1)' ;

period = 2*pi;
if stencil==0
    matrix = m20130226_06_FourierSpectralInterpolationMatrix(numel(theta1), theta2);
    if theta1(1) ~= 0
        error('The Fourier interpolation only works when the first element is 0.')
    end
else
    matrix = m20160925_01_periodicInterpolation(theta1,theta2,period,stencil);
end

figure(1)
clf

f1 = f(theta1);
plot([theta1;theta1+period],[f1;f1],'ob','DisplayName','Function on grid 1')
hold on
plot([theta2-period;theta2;theta2+period],[f(theta2);f(theta2);f(theta2)],'.-','Color',[0,0.7,0],'DisplayName','True function on grid 2')
temp = matrix*f1;
plot([theta2;theta2+period],[temp;temp],'.--r','DisplayName','Interpolated function on grid 2')
legend show
xlabel('theta')

assignin('base','matrix',matrix)

end