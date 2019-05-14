%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
function [qam_pulse] = qam_bh_pulse(pulse_length,carrier_freq,d_t) 

[bh_func, bh_diff] = bh_window(0:d_t:pulse_length,pulse_length);
qam_pulse = complex(cumtrapz(bh_diff.*sin(2*pi*carrier_freq*[0:d_t:pulse_length])),...
cumtrapz(bh_diff.*cos(2*pi*carrier_freq*[0:d_t:pulse_length])));
qam_pulse = pulse_length*d_t*qam_pulse;

end
