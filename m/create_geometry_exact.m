%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

pderect(0.4*[-1 1 -1 1],'SQ1'); 
pdecirc(0,0,0.16,'C1');
boundary_param = 64;
make_boundary_curve;
pdepoly(p(:,1),p(:,2),'P1');

n_phi = 20;
theta = pi/8 - pi/4;
phi = [2*pi/n_phi:2*pi/n_phi:2*pi];
r = (15+cos(2*phi+2)+sin(3*phi));  r = r/max(r);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
p =  repmat([0.07 ; 0.0],1,length(phi)) +0.018*R*[0.3*cos(phi); 2.5*sin(phi) ];
p = p';
pdepoly(p(:,1),p(:,2),'P2');

n_phi = 20;
theta = 7*pi/8 + pi;
phi = [2*pi/n_phi:2*pi/n_phi:2*pi];
r = (55+cos(5*phi+3)+sin(3*phi));  r = r/max(r);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
p =  repmat([-0.02 ; -0.06],1,length(phi)) +0.02*R*[1.7*cos(phi) + 0.2*cos(2*phi); 0.6*sin(phi) + 0.1*sin(2*phi+2)];
p = p';
pdepoly(p(:,1),p(:,2),'P3');

n_phi = 20;
theta = pi/2;
phi = [2*pi/n_phi:2*pi/n_phi:2*pi];
r = (15+cos(2*phi+5)+sin(1*phi));  r = r/max(r);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
p =  repmat([-0.04 ; 0.035],1,length(phi)) +0.018*R*[0.3*cos(phi) ; 2.4*sin(phi) ];
p = p';
pdepoly(p(:,1),p(:,2),'P4');

make_boundary_curve;
p = 0.92*p;
pdepoly(p(:,1),p(:,2),'P5');
