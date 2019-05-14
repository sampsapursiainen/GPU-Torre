%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
n_phi = 70;
phi = [2*pi/n_phi:2*pi/n_phi:2*pi];
r = (10+cos(5*phi+2)+sin(3*phi+1));    
p=[1.2*r.*cos(phi);1*r.*sin(phi)]';
p_norm = sqrt(sum(p.^2,2));
p = 0.145*p./max(p_norm);
c_x = cos([0 : 0.01 : 2*pi]);
c_y = sin([0 : 0.01 : 2*pi]);
