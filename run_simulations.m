%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
parameters;

for monobi = 1:2
n_constellation = monobi;

for tapaus=1:6
n_inv_iter = [1 2 3 1 2 3];
n_tv_iter = [3 3 3 3 3 3];
n_born = [1 1 1 2 2 2];

n_inv_iter = n_inv_iter(tapaus);
n_tv_iter = n_tv_iter(tapaus);
n_born = n_born(tapaus);

if monobi == 1
static = 'monostatic_';
else
static = 'bistatic_';
end
case_num = [static  num2str(tapaus)];

inversion_procedure;

['Done case ' case_num]
end
end
