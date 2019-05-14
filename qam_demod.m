%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
function [data_val,amp_val] = qam_demod(wave_val,carrier_freq,data_ind,data_param,d_t)

t_val = (data_ind - 1)*data_param*d_t;

carrier_1 = sin(2*pi*carrier_freq*t_val);
carrier_2 = cos(2*pi*carrier_freq*t_val);

data_val = real(wave_val).*carrier_1 + imag(wave_val).*carrier_2;
amp_val = sqrt(real(wave_val).^2 + imag(wave_val).^2);

end
