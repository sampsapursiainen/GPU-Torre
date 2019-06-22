%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
pderect(0.4*[-1 1 -1 1],'SQ1'); 
pdecirc(0,0,0.16,'C1');
boundary_param = 64;
make_boundary_curve;
pdepoly(p(:,1),p(:,2),'P1');
