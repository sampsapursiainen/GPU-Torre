%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
clear all
clc
load_field_data;

load data/mesh_1.mat;
load data/system_data_1.mat;

system_setting_index = 1;
parameters; 

[path_data] = create_constellation(n_r,n_constellation, bistatic_difference, mixing_rate, sparsity_factor, back_scattering, plot_constellation); 

u_data = u_data_background; 
f_data = f_data_background;
u_data_quad = u_data_background_quad; 
f_data_quad = f_data_background_quad;

make_jacobian;

n_ast = size(triangles_ast,1);
n_field = size(ast_p_ind,1);
n_path = size(path_data,1);
n_t = length(t_vec(1:data_param:end));

T_0_ind = find(t_vec(1:data_param:end) >= T_0, 1);
T_1_ind = find(t_vec(1:data_param:end) >= T_1, 1);

rec_data = zeros(n_path, 2*n_t);

for i = 1 : n_path
    rec_data(i,:) = rec_data_2(path_data(i,1),:,path_data(i,2)) - rec_data_1(path_data(i,1),:,path_data(i,2));  
end

TV_D = full(TV_matrix(nodes(ast_p_ind,:),triangles_ast));

x = zeros(n_ast,1);
theta = ones(size(x));
rng('default')

y = reshape(rec_data(:,T_0_ind:T_1_ind), [(T_1_ind-T_0_ind+1)*n_path 1]);
relative_noise = 20*log10(noise_level/max(abs(y)))
y = y + noise_level*randn(size(y));


for j = 1 : n_inv_iter

L = reshape(J_mat(:,T_0_ind:T_1_ind,:), [(T_1_ind-T_0_ind+1)*n_path n_ast]);

% Backpropagation
x_update = inv_omega*L'*y/norm(L,'fro');

y = y - L*x_update;
x = x + x_update;

plot_reconstruction;

if j < n_inv_iter
make_correction;
make_jacobian;
end

end
clear LTL;
error_evaluation