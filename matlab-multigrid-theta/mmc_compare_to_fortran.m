function mmc_compare_to_Fortran(filename)

global Ntheta Nzeta Nxi
global nu E
global theta_derivative_option preconditioner_theta_derivative_option
global zeta_derivative_option preconditioner_zeta_derivative_option
global xi_derivative_option preconditioner_xi_derivative_option
global pitch_angle_scattering_option preconditioner_pitch_angle_scattering_option
global xi_quadrature_option constraint_option coarsen_theta coarsen_zeta coarsen_xi
global Ntheta_min Nzeta_min Nxi_min smoothing_option restriction_option N_smoothing
global flux flow

fid = fopen(filename);
if fid<0
    fprintf('Unable to open %s to compare with fortran results.\n',filename)
    return
else
    fprintf('Comparing matlab results to fortran results in %s\n',filename)
end

% First compare 'input' quantities that should agree to within roundoff error.
tolerance = 1e-12;
comparisonType = 1;
% 1 = max(abs(difference)) < tolerance
% 2 = max(abs(difference) ./ mean) < tolerance

compare('Ntheta')
compare('Nzeta')
compare('Nxi')
compare('nu')
compare('E')
compare('theta_derivative_option')
compare('preconditioner_theta_derivative_option')
compare('zeta_derivative_option')
compare('preconditioner_zeta_derivative_option')
compare('xi_derivative_option')
compare('preconditioner_xi_derivative_option')
compare('pitch_angle_scattering_option')
compare('preconditioner_pitch_angle_scattering_option')
compare('xi_quadrature_option')
compare('constraint_option')
compare('coarsen_theta')
compare('coarsen_zeta')
compare('coarsen_xi')
compare('Ntheta_min')
compare('Nzeta_min')
compare('Nxi_min')
compare('smoothing_option')
compare('restriction_option')
compare('N_smoothing')
  
% Then compare 'output' quantities that only will agree within the
% convergence tolerance:
tolerance = 1e-3;
comparisonType = 2;
% 1 = max(abs(difference)) < tolerance
% 2 = max(abs(difference) ./ mean) < tolerance

compare('flux')
compare('flow')

fclose(fid);



    function compare(matlab_name)
        matlab_value_1 = eval(matlab_name);
        if islogical(matlab_value_1)
            if matlab_value_1
                matlab_value=1;
            else
                matlab_value=0;
            end
        else
            matlab_value = matlab_value_1;
        end
        
        line = fgetl(fid);
        [fortran_value, count] = sscanf(line,'%g');
        if count ~= 1
            fprintf('** Error parsing the following line from %s:\n',filename)
            fprintf([line,'\n'])
            return
        end
        difference = abs(matlab_value - fortran_value);
        switch comparisonType
            case 1
                difference_to_compare = difference;
            case 2
                avg = abs(matlab_value + fortran_value)/2;
                difference_to_compare = difference / avg;
            otherwise
                error('Invalid comparisonType')
        end
        
        if difference_to_compare < tolerance
            fprintf('  %s agrees. Matlab = %g, Fortran = %g\n',matlab_name,matlab_value,fortran_value)
        else
            fprintf('** ERROR!  %s disagrees! scalarDifference = %g\n',matlab_name, difference_to_compare)
            fprintf('  Matlab version: %g,  Fortran version: %g\n',matlab_value, fortran_value)
        end
    end

end