function [grad_vec] = grad_prod(B,u,triangles)

grad_vec = double(sum(B.*u(triangles),2));

end 


