%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

function [grad_vec] = grad_prod(B,u,triangles)

grad_vec = double(sum(B.*u(triangles),2));

end 


