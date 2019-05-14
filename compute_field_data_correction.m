%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
load data/mesh_1.mat
load data/system_data_1.mat

i_ind = 1;
parameters;

d_ind_aux_vec = ast_p_ind(:);
data_name = 'field_data_correction';
source_list = 1 : length(d_ind_aux_vec);
if exist('n_processes')
source_list = source_list(process_num : n_processes : end);
if exist('process_gpu_num')
gpu_num = process_gpu_num;
gpuDevice(gpu_num);
end
end

field_iteration;

