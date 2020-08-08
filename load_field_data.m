%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

load data/mesh_1.mat;
load data/system_data_1.mat;

i_ind = 1;
parameters;

n_p = size(ast_p_ind,1);
n_t_bg = length(t_vec(1:data_param:end));

u_data_background = zeros(n_r, n_t_bg, n_p);
f_data_background = zeros(n_r, n_t_bg, n_p); 
u_data_correction = zeros(n_r, n_t_bg, n_p);
f_data_correction = zeros(n_r, n_t_bg, n_p); 
u_data_background_quad = zeros(n_r, n_t_bg, n_p);
f_data_background_quad = zeros(n_r, n_t_bg, n_p); 
u_data_correction_quad = zeros(n_r, n_t_bg, n_p);
f_data_correction_quad = zeros(n_r, n_t_bg, n_p); 

rec_data_1 = zeros(n_r, 2*n_t_bg, n_r); 
rec_data_2 = zeros(n_r, 2*n_t_bg, n_r); 

data_ind_aux = 0;

 for k = 1 : n_r
 load(['field_data_background/u_data_' int2str(k) '.mat']); 
 load(['field_data_background/f_data_' int2str(k) '.mat']); 
 u_data_background(k, :, :) = u_data;
 f_data_background(k, :, :) = f_data;
 u_data_background_quad(k, :, :) = u_data_quad;
 f_data_background_quad(k, :, :) = f_data_quad;
 rec_data_1(k,:,:) = [rec_data;rec_data_quad];
 end

  for k = 1 : n_r
 load(['field_data_exact/u_data_' int2str(k) '.mat']); 
 rec_data_2(k,:,:) = [rec_data;rec_data_quad];
  end
 
 if n_inv_iter > 1

 for k = 1 : n_p
 load(['field_data_correction/u_data_' int2str(k) '.mat']); 
 load(['field_data_correction/f_data_' int2str(k) '.mat']); 
 u_data_correction(k, :, :) = u_data;
 f_data_correction(k, :, :) = f_data; 
 u_data_correction_quad(k, :, :) = u_data_quad;
 f_data_correction_quad(k, :, :) = f_data_quad; 
 end
 
 end
 

