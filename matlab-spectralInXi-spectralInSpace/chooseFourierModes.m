function [ms,ns,Binv_vec] = chooseFourierModes(geometryParameters,resolutionParameters)

NFourier = resolutionParameters.NFourier;

if geometryParameters.axisymmetric
    % Axisymmetric, so just include the first NFourier poloidal modes:
    ms = 0:(NFourier-1)';
    ns = zeros(size(ms));
else
    
    mmax = 32;
    nmax = 32;
    Ntheta = mmax*2+3;
    Nzeta = nmax*2+3;
    
    % Build theta grid:
    theta = linspace(0,2*pi,Ntheta+1);
    theta(end)=[];
    
    zetaMax = 2*pi/geometryParameters.Nperiods;
    zeta = linspace(0,zetaMax,Nzeta+1);
    zeta(end)=[];
    
    [zeta2D, theta2D] = meshgrid(zeta,theta);
    
    Nperiods = geometryParameters.Nperiods;
    epsilon_t = geometryParameters.epsilon_t;
    epsilon_h = geometryParameters.epsilon_h;
    helicity_l = geometryParameters.helicity_l;
    
    B = 1 + epsilon_t * cos(theta2D) + epsilon_h * cos(helicity_l*theta2D-Nperiods*zeta2D);
    
    % Add B to ensure all the modes of B have a shot at being included.
    BPower = B .^ (-2) + B;
    Binv = 1./B;
    
    % (slow) Fourier transform
    m1D = 0:mmax;
    n1D = (-nmax):nmax;
    [n2D, m2D] = meshgrid(n1D,m1D);
    amplitudes = zeros(size(n2D));
    Binv_Fourier2D = zeros(size(n2D));
    if NFourier>numel(n2D)
        error('Requested NFourier exceeds the number of modes in mmax and nmax.\n')
    end
    
    for m=0:mmax
        if m==0
            nmin=0;
        else
            nmin = -nmax;
        end
        for n=nmin:nmax
            if m==0 && n==0
                factor=1;
            else
                factor=2;
            end
            angle = m*theta2D-n*zeta2D*Nperiods;
            sinangle = sin(angle);
            cosangle = cos(angle);
            amplitude_sin = abs(sum(sum(sinangle.*BPower)));
            amplitude_cos = abs(sum(sum(cosangle.*BPower)));
            amplitudes(m+1,n+nmax+1) = amplitude_sin*amplitude_sin + amplitude_cos*amplitude_cos;
            Binv_Fourier2D(m+1,n+nmax+1) = factor*sum(sum(cosangle.*Binv));
        end
    end
    Binv_Fourier2D = Binv_Fourier2D / (Ntheta*Nzeta);
    
    m=0;n=0;amplitudes(m+1,n+nmax+1) = 0; % Zero out the amplitude for m=n=0, since we will manually put this mode in the first element of the ms and ns arrays.
    
    
    s = numel(m2D);
    m2D_reshaped = reshape(m2D,[s,1]);
    n2D_reshaped = reshape(n2D,[s,1]);
    amplitudes_reshaped = reshape(amplitudes,[s,1]);
    Binv_Fourier2D_reshaped = reshape(Binv_Fourier2D,[s,1]);
    [amplitudes_sorted, permutation] = sort(amplitudes_reshaped,'descend');
    m2D_sorted = m2D_reshaped(permutation);
    n2D_sorted = n2D_reshaped(permutation);
    ms = [0;m2D_sorted(1:NFourier-1)];
    ns = [0;n2D_sorted(1:NFourier-1)];
    fprintf('Smallest amplitude included: %g\n',amplitudes_sorted(NFourier-1))
    
    % Sort for increasing n:
    [ns,permutation2] = sort(ns);
    ms = ms(permutation2);

    % Now sort for increasing m:
    [ms,permutation3] = sort(ms);
    ns = ns(permutation3);
    
    Binv_vec = zeros(NFourier*2-1,1);
    %Binv_vec(1:NFourier) = Binv_Fourier2D(ms+1, ns+nmax+1);
    for imn = 1:NFourier
        index = find(ms(imn)==m2D_reshaped & ns(imn)==n2D_reshaped);
        if numel(index) ~= 1
            error('Should not get here!')
        else
            Binv_vec(imn) = Binv_Fourier2D_reshaped(index);
        end
    end
    
    included = zeros(size(m2D));
    for i=1:resolutionParameters.NFourier
        included(m2D==ms(i) & n2D==ns(i))=1;
    end
    
    figure(20)
    clf
    numRows=2;
    numCols=2;
    
    subplot(numRows,numCols,1)
    imagesc([1,Nzeta],[1,Ntheta],B)
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title('B')
    
    subplot(numRows,numCols,2)
    imagesc([1,Nzeta],[1,Ntheta],BPower)
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title('BPower')
    
    subplot(numRows,numCols,3)
    %imagesc([0,mmax],[-nmax,nmax],log10(amplitudes))
    imagesc(log10(amplitudes))
    colorbar
    xlabel('n')
    ylabel('m')
    title('log10(amplitudes)')
    
    subplot(numRows,numCols,4)
    imagesc([-nmax,nmax],[0,mmax],included)
    colorbar
    xlabel('n')
    ylabel('m')
    title('included')
    
    ns = ns * Nperiods;
    
end