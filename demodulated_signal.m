%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
function [demodulated_signal_vec] = demodulated_signal(modulated_signal_vec, carrier_freq, d_t)

modulated_signal_vec = modulated_signal_vec(:)';
freq_ind = round(carrier_freq*d_t*length(modulated_signal_vec-1));
signal_aux = fft(modulated_signal_vec);
norm_1 = sqrt(sum(abs(modulated_signal_vec-mean(modulated_signal_vec)).^2));
signal_aux = [0 signal_aux(freq_ind+2:1+round(length(signal_aux-1)/2))...
    zeros(1,2*freq_ind) signal_aux(round(length(signal_aux-1)/2)+2:end-freq_ind)];
demodulated_signal_vec = real(ifft(signal_aux));
norm_2 = sqrt(sum(abs(demodulated_signal_vec.^2)));
demodulated_signal_vec = norm_1*demodulated_signal_vec/norm_2;

demodulated_signal_vec = demodulated_signal_vec - demodulated_signal_vec(1);

end