%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre


n_inv_iter = 1;
n_tv_iter = 1;
n_born = 1;
inv_alpha = 0.2;
inv_beta = 0.001;

n_r = 16;
n_constellation = 1;
mixing_rate = 1;
sparsity_factor = 1;
back_scattering = 1;
plot_constellation = 1;
bistatic_difference = 1;

boundary_param = 640;

noise_level = 0.05;

T = 1.1;
T_0 = 0.1;
T_1 = 1.1;

d_reg_param = 1e-2;

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
