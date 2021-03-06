%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

system_setting_index = 1;
parameters;

n_ast = size(triangles_ast,1);
n_field = size(ast_p_ind,1);
n_path = size(path_data,1);
n_t = length(t_vec(1:data_param:end));

[bh_pulse,d_bh_pulse] = bh_window(t_vec(1:data_param:end),pulse_length,0,'complex');
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

   Aux_mat = [nodes_ast(:,triangles_ast(i,1)) ; nodes_ast(:,triangles_ast(i,2))] - repmat(nodes_ast(:,triangles_ast(i,3)),2,1); 
    ala = abs((Aux_mat(ind_m(1,1),:).*Aux_mat(ind_m(2,2),:)-Aux_mat(ind_m(1,2),:).*Aux_mat(ind_m(2,1),:)))/2;
    

if use_gpu
f_data_aux = gpuArray(f_data);
u_data_aux = gpuArray(u_data);
f_data_aux_correction = gpuArray(f_data_correction);
u_data_aux_correction = gpuArray(u_data_correction);
zeros_aux_1 = gpuArray(zeros(size(f_data_aux,1), size(f_data_aux,2)));
zeros_aux_2 = gpuArray(zeros(size(u_data_aux_correction,1), size(u_data_aux_correction,2)));
else
f_data_aux = f_data;
u_data_aux = u_data;
f_data_aux_correction = f_data_correction;
u_data_aux_correction = u_data_correction;
zeros_aux_1 = zeros(size(f_data_aux,1), size(f_data_aux,2));
zeros_aux_2 = zeros(size(u_data_aux_correction,1), size(u_data_aux_correction,2));
end

nodes_ast = nodes(ast_p_ind,:)';

for i = 1 : n_ast
    for k = 1 : n_born
    
    fft_aux_f_path_1_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,1)) zeros_aux_1], [], 2);
    fft_aux_f_path_2_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,2)) zeros_aux_1], [], 2);
    fft_aux_f_path_3_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,3)) zeros_aux_1], [], 2);
    
    fft_aux_f_path_1 = ala*(D_mat(1,1)*fft_aux_f_path_1_aux + D_mat(1,2)*fft_aux_f_path_2_aux +  D_mat(1,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_2 = ala*(D_mat(2,1)*fft_aux_f_path_1_aux + D_mat(2,2)*fft_aux_f_path_2_aux +  D_mat(2,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_3 = ala*(D_mat(3,1)*fft_aux_f_path_1_aux + D_mat(3,2)*fft_aux_f_path_2_aux +  D_mat(3,3)*fft_aux_f_path_3_aux); 
    
    if correct_full_field
    
    fft_aux_f_path_1_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,1)) zeros_aux_2], [], 2);
    fft_aux_f_path_2_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,2)) zeros_aux_2], [], 2);
    fft_aux_f_path_3_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,3)) zeros_aux_2], [], 2);
    
    fft_aux_f_path_1_correction = ala*(D_mat(1,1)*fft_aux_f_path_1_aux + D_mat(1,2)*fft_aux_f_path_2_aux +  D_mat(1,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_2_correction = ala*(D_mat(2,1)*fft_aux_f_path_1_aux + D_mat(2,2)*fft_aux_f_path_2_aux +  D_mat(2,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_3_correction = ala*(D_mat(3,1)*fft_aux_f_path_1_aux + D_mat(3,2)*fft_aux_f_path_2_aux +  D_mat(3,3)*fft_aux_f_path_3_aux);
    
    end
    
    fft_aux_u_field_1 = fft([zeros_aux_2 u_data_aux_correction(:, :, triangles_ast(i,1)) zeros_aux_2],[],2);
    fft_aux_u_field_2 = fft([zeros_aux_2 u_data_aux_correction(:, :, triangles_ast(i,2)) zeros_aux_2],[],2);
    fft_aux_u_field_3 = fft([zeros_aux_2 u_data_aux_correction(:, :, triangles_ast(i,3)) zeros_aux_2],[],2);
    
    fft_aux_f_field_1_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,1)) zeros_aux_2],[],2);
    fft_aux_f_field_2_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,2)) zeros_aux_2],[],2);
    fft_aux_f_field_3_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,3)) zeros_aux_2],[],2);

    fft_aux_f_field_1 = ala*(D_mat(1,1)*fft_aux_f_field_1_aux + D_mat(1,2)*fft_aux_f_field_2_aux +  D_mat(1,3)*fft_aux_f_field_3_aux);
    fft_aux_f_field_2 = ala*(D_mat(2,1)*fft_aux_f_field_1_aux + D_mat(2,2)*fft_aux_f_field_2_aux +  D_mat(2,3)*fft_aux_f_field_3_aux);
    fft_aux_f_field_3 = ala*(D_mat(3,1)*fft_aux_f_field_1_aux + D_mat(3,2)*fft_aux_f_field_2_aux +  D_mat(3,3)*fft_aux_f_field_3_aux);
    
deconv_vec_u_1 = fft_aux_u_field_1./(fft_aux+d_reg_param);
deconv_vec_u_2 = fft_aux_u_field_2./(fft_aux+d_reg_param);
deconv_vec_u_3 = fft_aux_u_field_3./(fft_aux+d_reg_param);

deconv_vec_f_1 = fft_aux_f_field_1./(fft_aux+d_reg_param);
deconv_vec_f_2 = fft_aux_f_field_2./(fft_aux+d_reg_param);
deconv_vec_f_3 = fft_aux_f_field_3./(fft_aux+d_reg_param);
    
 for j = 1 : n_r
  
  aux_vec_1 = real(ifft(deconv_vec_u_1.*fft_aux_f_path_1(j,:),[],2));
  aux_vec_2 = real(ifft(deconv_vec_u_2.*fft_aux_f_path_2(j,:),[],2));
  aux_vec_3 = real(ifft(deconv_vec_u_3.*fft_aux_f_path_3(j,:),[],2));  
  aux_data = aux_vec_1(:,n_t+1:2*n_t) + aux_vec_2(:,n_t+1:2*n_t) + aux_vec_3(:,n_t+1:2*n_t);
  u_data_aux(j,:,:) = u_data_aux(j,:,:) + reshape(x_update(i)*aux_data',[1 n_t n_field]); 
    
  aux_vec_1 = real(ifft(deconv_vec_f_1.*fft_aux_f_path_1(j,:),[],2));
  aux_vec_2 = real(ifft(deconv_vec_f_2.*fft_aux_f_path_2(j,:),[],2));
  aux_vec_3 = real(ifft(deconv_vec_f_3.*fft_aux_f_path_3(j,:),[],2));  
  aux_data = aux_vec_1(:,n_t+1:2*n_t) + aux_vec_2(:,n_t+1:2*n_t) + aux_vec_3(:,n_t+1:2*n_t); 
  f_data_aux(j,:,:) = f_data_aux(j,:,:) + reshape(x_update(i)*aux_data',[1 n_t n_field]); 

 end
 
 if correct_full_field

  for j = 1 : n_p
  
  aux_vec_1 = real(ifft(deconv_vec_u_1.*fft_aux_f_path_1_correction(j,:),[],2));
  aux_vec_2 = real(ifft(deconv_vec_u_2.*fft_aux_f_path_2_correction(j,:),[],2));
  aux_vec_3 = real(ifft(deconv_vec_u_3.*fft_aux_f_path_3_correction(j,:),[],2));  
  aux_data = aux_vec_1(:,n_t+1:2*n_t) + aux_vec_2(:,n_t+1:2*n_t) + aux_vec_3(:,n_t+1:2*n_t);
  u_data_aux_correction(j,:,:) = u_data_aux_correction(j,:,:) + reshape(x_update(i)*aux_data',[1 n_t n_field]); 
    
  aux_vec_1 = real(ifft(deconv_vec_f_1.*fft_aux_f_path_1_correction(j,:),[],2));
  aux_vec_2 = real(ifft(deconv_vec_f_2.*fft_aux_f_path_2_correction(j,:),[],2));
  aux_vec_3 = real(ifft(deconv_vec_f_3.*fft_aux_f_path_3_correction(j,:),[],2));  
  aux_data = aux_vec_1(:,n_t+1:2*n_t) + aux_vec_2(:,n_t+1:2*n_t) + aux_vec_3(:,n_t+1:2*n_t); 
  f_data_aux_correction(j,:,:) = f_data_aux_correction(j,:,:) + reshape(x_update(i)*aux_data',[1 n_t n_field]); 

  end

 end
  
  time_val = toc;
  if mod(i,25)==0
   waitbar(i/(2*n_ast),h,['Deconvolution. Ready approx: ' datestr(datevec(now+(2*n_ast/i - 1)*time_val/86400)) '.']);
  end
    end
end

if use_gpu
f_data = gather(f_data_aux);
u_data = gather(u_data_aux);
f_data_correction = gather(f_data_aux_correction);
u_data_correction = gather(u_data_aux_correction);
else
f_data = f_data_aux;
u_data = u_data_aux;
f_data_correction = f_data_aux_correction;
u_data_correction = u_data_aux_correction;
end


if use_gpu
f_data_aux = gpuArray(f_data_quad);
u_data_aux = gpuArray(u_data_quad);
f_data_aux_correction = gpuArray(f_data_correction_quad);
u_data_aux_correction = gpuArray(u_data_correction_quad);
zeros_aux_1 = gpuArray(zeros(size(f_data_aux,1), size(f_data_aux,2)));
zeros_aux_2 = gpuArray(zeros(size(u_data_aux_correction,1), size(u_data_aux_correction,2)));
else
f_data_aux = f_data_quad;
u_data_aux = u_data_quad;
f_data_aux_correction = f_data_correction_quad;
u_data_aux_correction = u_data_correction_quad;
zeros_aux_1 = zeros(size(f_data_aux,1), size(f_data_aux,2));
zeros_aux_2 = zeros(size(u_data_aux_correction,1), size(u_data_aux_correction,2));
end

nodes_ast = nodes(ast_p_ind,:)';

for i = 1 : n_ast
    for k = 1 : n_born
    
    fft_aux_f_path_1_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,1)) zeros_aux_1], [], 2);
    fft_aux_f_path_2_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,2)) zeros_aux_1], [], 2);
    fft_aux_f_path_3_aux = fft([zeros_aux_1 f_data_aux(:, :, triangles_ast(i,3)) zeros_aux_1], [], 2);
    
    fft_aux_f_path_1 = ala*(D_mat(1,1)*fft_aux_f_path_1_aux + D_mat(1,2)*fft_aux_f_path_2_aux +  D_mat(1,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_2 = ala*(D_mat(2,1)*fft_aux_f_path_1_aux + D_mat(2,2)*fft_aux_f_path_2_aux +  D_mat(2,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_3 = ala*(D_mat(3,1)*fft_aux_f_path_1_aux + D_mat(3,2)*fft_aux_f_path_2_aux +  D_mat(3,3)*fft_aux_f_path_3_aux); 
    
    if correct_full_field
    
    fft_aux_f_path_1_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,1)) zeros_aux_2], [], 2);
    fft_aux_f_path_2_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,2)) zeros_aux_2], [], 2);
    fft_aux_f_path_3_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,3)) zeros_aux_2], [], 2);
    
    fft_aux_f_path_1_correction = ala*(D_mat(1,1)*fft_aux_f_path_1_aux + D_mat(1,2)*fft_aux_f_path_2_aux +  D_mat(1,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_2_correction = ala*(D_mat(2,1)*fft_aux_f_path_1_aux + D_mat(2,2)*fft_aux_f_path_2_aux +  D_mat(2,3)*fft_aux_f_path_3_aux);
    fft_aux_f_path_3_correction = ala*(D_mat(3,1)*fft_aux_f_path_1_aux + D_mat(3,2)*fft_aux_f_path_2_aux +  D_mat(3,3)*fft_aux_f_path_3_aux);
    
    end
    
    fft_aux_u_field_1 = fft([zeros_aux_2 u_data_aux_correction(:, :, triangles_ast(i,1)) zeros_aux_2],[],2);
    fft_aux_u_field_2 = fft([zeros_aux_2 u_data_aux_correction(:, :, triangles_ast(i,2)) zeros_aux_2],[],2);
    fft_aux_u_field_3 = fft([zeros_aux_2 u_data_aux_correction(:, :, triangles_ast(i,3)) zeros_aux_2],[],2);
    
    fft_aux_f_field_1_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,1)) zeros_aux_2],[],2);
    fft_aux_f_field_2_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,2)) zeros_aux_2],[],2);
    fft_aux_f_field_3_aux = fft([zeros_aux_2 f_data_aux_correction(:, :, triangles_ast(i,3)) zeros_aux_2],[],2);

    fft_aux_f_field_1 = ala*(D_mat(1,1)*fft_aux_f_field_1_aux + D_mat(1,2)*fft_aux_f_field_2_aux +  D_mat(1,3)*fft_aux_f_field_3_aux);
    fft_aux_f_field_2 = ala*(D_mat(2,1)*fft_aux_f_field_1_aux + D_mat(2,2)*fft_aux_f_field_2_aux +  D_mat(2,3)*fft_aux_f_field_3_aux);
    fft_aux_f_field_3 = ala*(D_mat(3,1)*fft_aux_f_field_1_aux + D_mat(3,2)*fft_aux_f_field_2_aux +  D_mat(3,3)*fft_aux_f_field_3_aux);
    
deconv_vec_u_1 = fft_aux_u_field_1./(fft_aux+d_reg_param);
deconv_vec_u_2 = fft_aux_u_field_2./(fft_aux+d_reg_param);
deconv_vec_u_3 = fft_aux_u_field_3./(fft_aux+d_reg_param);

deconv_vec_f_1 = fft_aux_f_field_1./(fft_aux+d_reg_param);
deconv_vec_f_2 = fft_aux_f_field_2./(fft_aux+d_reg_param);
deconv_vec_f_3 = fft_aux_f_field_3./(fft_aux+d_reg_param);
    
 
 for j = 1 : n_r
  
  aux_vec_1 = real(ifft(deconv_vec_u_1.*fft_aux_f_path_1(j,:),[],2));
  aux_vec_2 = real(ifft(deconv_vec_u_2.*fft_aux_f_path_2(j,:),[],2));
  aux_vec_3 = real(ifft(deconv_vec_u_3.*fft_aux_f_path_3(j,:),[],2));  
  aux_data = aux_vec_1(:,n_t+1:2*n_t) + aux_vec_2(:,n_t+1:2*n_t) + aux_vec_3(:,n_t+1:2*n_t);
  u_data_aux(j,:,:) = u_data_aux(j,:,:) + reshape(x_update(i)*aux_data',[1 n_t n_field]); 
    
  aux_vec_1 = real(ifft(deconv_vec_f_1.*fft_aux_f_path_1(j,:),[],2));
  aux_vec_2 = real(ifft(deconv_vec_f_2.*fft_aux_f_path_2(j,:),[],2));
  aux_vec_3 = real(ifft(deconv_vec_f_3.*fft_aux_f_path_3(j,:),[],2));  
  aux_data = aux_vec_1(:,n_t+1:2*n_t) + aux_vec_2(:,n_t+1:2*n_t) + aux_vec_3(:,n_t+1:2*n_t); 
  f_data_aux(j,:,:) = f_data_aux(j,:,:) + reshape(x_update(i)*aux_data',[1 n_t n_field]); 

 end

 if correct_full_field
 
  for j = 1 : n_p
  
  aux_vec_1 = real(ifft(deconv_vec_u_1.*fft_aux_f_path_1_correction(j,:),[],2));
  aux_vec_2 = real(ifft(deconv_vec_u_2.*fft_aux_f_path_2_correction(j,:),[],2));
  aux_vec_3 = real(ifft(deconv_vec_u_3.*fft_aux_f_path_3_correction(j,:),[],2));  
  aux_data = aux_vec_1(:,n_t+1:2*n_t) + aux_vec_2(:,n_t+1:2*n_t) + aux_vec_3(:,n_t+1:2*n_t);
  u_data_aux_correction(j,:,:) = u_data_aux_correction(j,:,:) + reshape(x_update(i)*aux_data',[1 n_t n_field]); 
    
  aux_vec_1 = real(ifft(deconv_vec_f_1.*fft_aux_f_path_1_correction(j,:),[],2));
  aux_vec_2 = real(ifft(deconv_vec_f_2.*fft_aux_f_path_2_correction(j,:),[],2));
  aux_vec_3 = real(ifft(deconv_vec_f_3.*fft_aux_f_path_3_correction(j,:),[],2));  
  aux_data = aux_vec_1(:,n_t+1:2*n_t) + aux_vec_2(:,n_t+1:2*n_t) + aux_vec_3(:,n_t+1:2*n_t); 
  f_data_aux_correction(j,:,:) = f_data_aux_correction(j,:,:) + reshape(x_update(i)*aux_data',[1 n_t n_field]); 

  end

 end
  
  time_val = toc;
  if mod(i,25)==0
   waitbar((n_ast+i)/(2*n_ast),h,['Deconvolution. Ready approx: ' datestr(datevec(now+(2*n_ast/(i+n_ast) - 1)*time_val/86400)) '.']);
  end
    end
end

if use_gpu
f_data_quad = gather(f_data_aux);
u_data_quad = gather(u_data_aux);
f_data_correction_quad = gather(f_data_aux_correction);
u_data_correction_quad = gather(u_data_aux_correction);
else
f_data_quad = f_data_aux;
u_data_quad = u_data_aux;
f_data_correction_quad = f_data_aux_correction;
u_data_correction_quad = u_data_aux_correction;
end

close(h);










