function mmc_defaults()

global epsilon_t epsilon_h iota G I Nperiods helicity_l
epsilon_t = -0.07053;
epsilon_h = 0.05067;
iota = 0.4542;
G = 3.7481;
I = 0;
Nperiods = 10;
helicity_l = 2;

global nu E
nu = 0.01;
E = 0;

global Ntheta Nzeta Nxi
Ntheta = 15;
Nzeta = 15;
Nxi = 30;

global Ntheta_min Nzeta_min Nxi_min
Ntheta_min = 7;
Nzeta_min = 7;
Nxi_min = 9;

global coarsen_theta coarsen_zeta coarsen_xi coarsen_option
coarsen_theta = false;
coarsen_zeta  = true;
coarsen_xi    = false;
coarsen_option = 1;

global theta_derivative_option zeta_derivative_option xi_derivative_option pitch_angle_scattering_option xi_quadrature_option
theta_derivative_option = 8;
zeta_derivative_option = 8;
xi_derivative_option = 8;
pitch_angle_scattering_option = 3;
xi_quadrature_option = 3;

global preconditioner_theta_derivative_option preconditioner_zeta_derivative_option preconditioner_xi_derivative_option preconditioner_pitch_angle_scattering_option
preconditioner_theta_derivative_option = 4;
preconditioner_zeta_derivative_option = 4;
preconditioner_xi_derivative_option = 4;
preconditioner_pitch_angle_scattering_option = 2;

global restriction_option
restriction_option = 1;
% 1 = all row sums equal 1 exactly
% 2 = restriction == prolongation' up to a constant.

global constraint_option
constraint_option = 1;

global Jacobi_omega N_smoothing smoothing_option
Jacobi_omega = 2/3;
N_smoothing = 1;
smoothing_option = 1;

global gmres_restart gmres_maxit gmres_tol
gmres_restart = 100;
gmres_maxit = 200;
gmres_tol = 1e-6;


end