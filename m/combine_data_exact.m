%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

system_setting_index = 2; 
parameters;

data_name = 'field_data_exact';

for k = 1 : n_r

k
load(['./' data_name '/u_data_modulated_' int2str(k) '.mat']);

[rec_data, rec_data_quad, rec_amp] = qam_demod(transpose(rec_data),carrier_freq,pulse_length,t_data);

rec_data = transpose(rec_data);
rec_data_quad = transpose(rec_data_quad);

save(['./' data_name '/u_data_' int2str(k) '.mat'], 'rec_data','rec_data_quad','t_data');

end

