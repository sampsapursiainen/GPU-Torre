function [div_vec] = div_prod(B_1,B_2,p_1,p_2,triangles)

aux_vec = B_1.*p_1 + B_2.*p_2;
div_vec = double(accumarray(triangles(:),aux_vec(:)));

end
