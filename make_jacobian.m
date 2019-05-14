%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
i_ind = 1;
parameters;

n_ast = size(triangles_ast,1);
n_field = size(ast_p_ind,1);
n_path = size(path_data,1);
n_t = length(t_vec(1:data_param:end));

[bh_pulse,d_bh_pulse] = bh_window([0:d_t:pulse_length],pulse_length);
signal_vec = zeros(1,length(t_vec));
signal_vec(1:length(bh_pulse)) = bh_pulse;
signal_vec = signal_vec(1:data_param:length(t_vec));
if use_gpu 
signal_vec = gpuArray([zeros(size(signal_vec)) signal_vec zeros(size(signal_vec))]);
else
signal_vec = [zeros(size(signal_vec)) signal_vec zeros(size(signal_vec))];
end

h = waitbar(0,'Deconvolution.');
tic;
fft_aux = fft(signal_vec,[],2);
D_mat = [2 1 1; 1 2 1; 1 1 2]/12;
ind_m = [1 3; 2 4];

if use_gpu
f_data_aux = gpuArray(f_data);
u_data_aux = gpuArray(u_data);
J_mat = gpuArray(zeros(n_path, n_t, n_ast));
zeros_aux_1 = gpuArray(zeros(size(f_data_aux,1), size(f_data_aux,2)));
else
f_data_aux = f_data;
u_data_aux = u_data;
J_mat = zeros(n_path, n_t, n_ast);
zeros_aux_1 = zeros(size(f_data_aux,1), size(f_data_aux,2));
end

nodes_ast = nodes(ast_p_ind,:)';

for i = 1 : n_ast
    for k = 1 : n_born
    
    Aux_mat = [nodes_ast(:,triangles_ast(i,1)) ; nodes_ast(:,triangles_ast(i,2))] - repmat(nodes_ast(:,triangles_ast(i,3)),2,1); 
    ala = abs((Aux_mat(ind_m(1,1),:).*Aux_mat(ind_m(2,2),:)-Aux_mat(ind_m(1,2),:).*Aux_mat(ind_m(2,1),:)))/2;
    
    fft_aux_f_path_1_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,1)) zeros_aux_1], [], 2);
    fft_aux_f_path_2_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,2)) zeros_aux_1], [], 2);
    fft_aux_f_path_3_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,3)) zeros_aux_1], [], 2);
    
    fft_aux_f_path_1 = ala*(D_mat(1,1)*fft_aux_f_path_1_aux + D_mat(1,2)*fft_aux_f_path_2_aux +  D_mat(1,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_2 = ala*(D_mat(2,1)*fft_aux_f_path_1_aux + D_mat(2,2)*fft_aux_f_path_2_aux +  D_mat(2,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_3 = ala*(D_mat(3,1)*fft_aux_f_path_1_aux + D_mat(3,2)*fft_aux_f_path_2_aux +  D_mat(3,3)*fft_aux_f_path_3_aux);

    fft_aux_u_path_1 = fft([zeros_aux_1 u_data_aux(:, :, triangles_ast(i,1)) zeros_aux_1], [], 2);
    fft_aux_u_path_2 = fft([zeros_aux_1 u_data_aux(:, :, triangles_ast(i,2)) zeros_aux_1], [], 2);
    fft_aux_u_path_3 = fft([zeros_aux_1 u_data_aux(:, :, triangles_ast(i,3)) zeros_aux_1], [], 2);
  
  for  j = 1 : n_path
    
    deconv_vec_1 = fft_aux_u_path_1(path_data(j,1),:)./(fft_aux+d_reg_param);
    deconv_vec_2 = fft_aux_u_path_2(path_data(j,1),:)./(fft_aux+d_reg_param);
    deconv_vec_3 = fft_aux_u_path_3(path_data(j,1),:)./(fft_aux+d_reg_param);
    aux_data_vec_1 = fft_aux_f_path_1(path_data(j,2),:);
    aux_data_vec_2 = fft_aux_f_path_2(path_data(j,2),:);
    aux_data_vec_3 = fft_aux_f_path_3(path_data(j,2),:);
          
    aux_vec_1 = real(ifft(deconv_vec_1.*aux_data_vec_1,[],2));
    aux_vec_2 = real(ifft(deconv_vec_2.*aux_data_vec_2,[],2));
    aux_vec_3 = real(ifft(deconv_vec_3.*aux_data_vec_3,[],2));
    aux_data = sum(aux_vec_1(:,n_t + 1 : 2*n_t) + aux_vec_2(:,n_t + 1 : 2*n_t) + aux_vec_3(:,n_t + 1 : 2*n_t),1);
    J_mat(j, :, i) = J_mat(j, :, i) + aux_data; 
    
  end

  
  time_val = toc;
  if mod(i,25)==0
   waitbar(i/size(triangles_ast,1),h,['Deconvolution. Ready approx: ' datestr(datevec(now+(size(triangles_ast,1)/i - 1)*time_val/86400)) '.']);
  end
    end
end

if use_gpu
J_mat = gather(J_mat);
end

close(h);










