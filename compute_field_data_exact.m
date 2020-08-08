%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

load data/mesh_2.mat
load data/system_data_2.mat

i_ind = 2;
parameters;

d_ind_aux_vec = r_ind(:);
data_name = 'field_data_exact';
source_list = 1 : length(d_ind_aux_vec);
if exist('n_processes')
source_list = source_list(process_num : n_processes : end);
if exist('process_gpu_num')
gpu_num = process_gpu_num;
gpuDevice(gpu_num);
end
end

field_iteration;
combine_data_exact;

