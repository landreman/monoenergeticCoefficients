function [B, dBdtheta, dBdzeta] = mmc_geometry(theta,zeta)

global epsilon_t epsilon_h Nperiods helicity_l

B = 1 + epsilon_t * cos(theta) + epsilon_h * cos(helicity_l*theta-Nperiods*zeta);

dBdtheta = -epsilon_t * sin(theta) - helicity_l * epsilon_h * sin(helicity_l*theta-Nperiods*zeta);

dBdzeta = Nperiods * epsilon_h * sin(helicity_l*theta-Nperiods*zeta);

end