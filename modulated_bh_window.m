%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
function [signal_vec] = modulated_bh_window(pulse_length, carrier_freq, d_t)

[bh_pulse] = bh_window(0 : d_t : pulse_length,pulse_length);

carrier_wave = zeros(size(bh_pulse));
freq_ind = round(carrier_freq*d_t*length(bh_pulse-1));
if freq_ind > 0 
carrier_wave(1+freq_ind) = 1;
carrier_wave(end-freq_ind+1) = 1;
else 
carrier_wave(1) = 2;   
end

carrier_wave = real(ifft(carrier_wave)*length(carrier_wave)/2);

signal_vec = bh_pulse.*carrier_wave;

end