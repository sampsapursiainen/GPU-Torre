%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

n_inv_iter = 1;
n_tv_iter = 10;
n_tkh_iter = 1;
n_born = 1;
inv_alpha = 0.5;
inv_beta = 0.005;%0.01
inv_omega = 0.01;
inv_epsilon = 1E-15;%0.01

correct_full_field = 0;

noise_level = 0.1; %noise level added to the data

start_point=1E-4;  %gradient descent iteration for inv_beta start point
end_point=1E4;     %gradient descent iteration  for inv_beta end point
discrepancy_h = 0.1; % discrepancy error for discretisation
latent_noise_factor = 1; %latent noise related to discretisation
N=50;                   % number of gradient descent sample point
n_lin_samp=10;          % number of linear sampling iteration
n_bisection=5;          % number of gradient descent iteration
    
n_r = 32;
n_constellation = 1;
mixing_rate = 1;
sparsity_factor = 1;
back_scattering = 1;
plot_constellation = 1;
bistatic_difference = 2;

plot_threshold_db = -30;

boundary_param = 640;


T = 1.1;
T_0 = 0.15;
T_1 = 0.85;

d_reg_param = 1e-4;

pml_val = 15;
or_radius = 0.16;

frame_param = 100;

ichol_tol = 1e-7;

use_gpu = 0;

pcg_tol = 1e-5;
pcg_max_it = 10000;

if i_ind == 1

pulse_length = 0.1;
carrier_freq = 0;   
d_t = 0.00025;
data_param = 20;
ast_aux_ind = [3];
bg_permittivity = 4;
absorb_val = 5*15.5/max(carrier_freq,15.5);
cav_ind = [];
cav_val = [];
n_ref_aux = 2;
gpu_num = 1;

elseif i_ind == 2
   
pulse_length = 0.1;  
carrier_freq = 0;    
d_t = 0.00025/4; 
data_param = 20*4;
ast_aux_ind = [2 3 5 6 7];
cover_aux_ind = [3]; 
cav_ind = [5 6 7];
cav_val = [1 1 1];
bg_permittivity = 4;
sl_permittivity = 3;
absorb_val = 5*15.5/max(carrier_freq,15.5);
n_ref_aux = 2;    
gpu_num = 2;

end

im_unit = complex(0,1);
carrier_omega = 2*pi*carrier_freq;
t_vec = [0:d_t:T];

if use_gpu 
gpuDevice(gpu_num);
end
