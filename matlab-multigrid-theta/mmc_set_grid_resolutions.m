global Ntheta Nzeta Nxi
global Nthetas Nzetas Nxis
global Ntheta_min Nzeta_min Nxi_min
global coarsen_theta coarsen_zeta coarsen_xi coarsen_option

if mod(Ntheta,2)==0
    Ntheta = Ntheta + 1;
end
if mod(Nzeta,2)==0
    Nzeta = Nzeta + 1;
end


small = 1.0e-12;
N_levels = 1;
if coarsen_theta
    N_levels = max([N_levels, (round(log(Ntheta/Ntheta_min) / log(2.0) - small) +1)]);
end
if coarsen_zeta
    N_levels = max([N_levels, (round(log(Nzeta/Nzeta_min) / log(2.0) - small) +1)]);
end
if coarsen_xi
    N_levels = max([N_levels, (round(log(Nxi/Nxi_min) / log(2.0) - small) +1)]);
end

Nthetas = zeros(N_levels,1);
Nzetas  = zeros(N_levels,1);
Nxis    = zeros(N_levels,1);

if N_levels==1
    Nthetas = Ntheta;
    Nzetas = Nzeta;
    Nxis = Nxi;
else
    if coarsen_theta
        for j=1:N_levels
            switch coarsen_option
                case 1
                    temp_float = max([Ntheta_min, Ntheta * (0.5 ^ (j-1))]);
                case 2
                    temp_float = exp(log(Ntheta) - (j-1)/(N_levels-1)*log(Ntheta/Ntheta_min));
                otherwise
                    error('Invalid coarsen_option')
            end
            
            temp_int = round(temp_float);
            if (mod(temp_int,2)==1)
                Nthetas(j) = temp_int;
            else
                % Check whether temp_int+1 or temp_int-1 is closer to the ideal value:
                if (abs(temp_int+1 - temp_float) > abs(temp_int-1 - temp_float))
                    Nthetas(j) = temp_int-1;
                else
                    Nthetas(j) = temp_int+1;
                end
            end
        end
    else
        Nthetas = Nthetas + Ntheta;
    end
    
    if coarsen_zeta
        for j=1:N_levels
            switch coarsen_option
                case 1
                    temp_float = max([Nzeta_min, Nzeta * (0.5 ^ (j-1))]);
                case 2
                    temp_float = exp(log(Nzeta) - (j-1)/(N_levels-1)*log(Nzeta/Nzeta_min));
                otherwise
                    error('Invalid coarsen_option')
            end
            
            temp_int = round(temp_float);
            if (mod(temp_int,2)==1)
                Nzetas(j) = temp_int;
            else
                % Check whether temp_int+1 or temp_int-1 is closer to the ideal value:
                if (abs(temp_int+1 - temp_float) > abs(temp_int-1 - temp_float))
                    Nzetas(j) = temp_int-1;
                else
                    Nzetas(j) = temp_int+1;
                end
            end
        end
    else
        Nzetas  = Nzetas  + Nzeta;
    end
    
    if coarsen_xi
        switch coarsen_option
            case 1
                for j=1:N_levels
                    Nxis(j) = max([Nxi_min, round((Nxi-1) * (0.5 ^ (j-1)) + 1)]);
                end
            case 2
                for j=1:N_levels
                    Nxis(j) = round(exp(log(Nxi-1) - (j-1)/(N_levels-1)*log((Nxi-1)/(Nxi_min-1)))) + 1;
                end
            otherwise
                error('Invalid coarsen_option')
        end
    else
        Nxis    = Nxis    + Nxi;
    end
end

fprintf('----- Computed parameters for multigrid: ----\n')
fprintf('  Level  Ntheta   Nzeta     Nxi\n')
for level = 1:N_levels
    fprintf('%7d %7d %7d %7d\n',level, Nthetas(level), Nzetas(level), Nxis(level))
end
