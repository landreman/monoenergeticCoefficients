function [B, dBdtheta, dBdzeta] = geometry(theta,zeta,geometryParameters)

Nperiods = geometryParameters.Nperiods;
epsilon_t = geometryParameters.epsilon_t;
epsilon_h = geometryParameters.epsilon_h;
helicity_l = geometryParameters.helicity_l;

B = 1 + epsilon_t * cos(theta) + epsilon_h * cos(helicity_l*theta-Nperiods*zeta);

dBdtheta = -epsilon_t * sin(theta) - helicity_l * epsilon_h * sin(helicity_l*theta-Nperiods*zeta);

dBdzeta = Nperiods * epsilon_h * sin(helicity_l*theta-Nperiods*zeta);

end