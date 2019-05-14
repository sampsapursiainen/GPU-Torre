%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
function [ bh_vec, d_bh_vec, d2_bh_vec ] = bh_window(t, T)

a_0 = 0.35875;
a_1 = 0.48829;
a_2 = 0.14128; 
a_3 = 0.01168;

bh_vec = a_0 - a_1*cos(2*pi*t/T) + a_2*cos(4*pi*t/T) - a_3*cos(6*pi*t/T);
d_bh_vec =  a_1*(2*pi/T)*sin(2*pi*t/T) - a_2*(4*pi/T)*sin(4*pi*t/T) + a_3*(6*pi/T)*sin(6*pi*t/T);
d2_bh_vec = a_1*(2*pi/T).^2*cos(2*pi*t/T) - a_2*(4*pi/T).^2*cos(4*pi*t/T) + a_3*(6*pi/T).^2*cos(6*pi*t/T);

a_mat = ones(size(t));
a_mat(t>T) = 0;
a_mat(t<0) = 0;
bh_vec = a_mat.*bh_vec;
d_bh_vec = a_mat.*d_bh_vec;
d2_bh_vec = a_mat.*d2_bh_vec;

end

