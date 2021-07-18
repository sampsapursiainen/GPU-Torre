%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

if system_setting_index == 1
n_p = length(ast_p_ind(:));
end

if not(use_gpu)
p_vec = symamd(C);
P_mat = speye(length(p_vec),length(p_vec));
i_p_vec = P_mat(:,p_vec)*[1:length(p_vec)]';
L_mat = ichol(C(p_vec,p_vec),struct('type','ict','droptol',ichol_tol,'michol','off'));
U_mat = L_mat';
end

B_1_T = B_1';
B_2_T = B_2';

if use_gpu
B_1 = gpuArray(B_1);
B_2 = gpuArray(B_2); 
B_1_T = gpuArray(B_1_T);
B_2_T = gpuArray(B_2_T); 
C = gpuArray(C);
R = gpuArray(R);
A = gpuArray((1./full(diag(A)))); 
M = gpuArray((1./full(sum(C,2))));
I_p_x = gpuArray(I_p_x);
I_p_y = gpuArray(I_p_y);
I_t_x = gpuArray(I_t_x);
I_t_y = gpuArray(I_t_y);
pml_val = gpuArray(pml_val);
end

make_boundary_curve;

figure(1); clf;
set(1,'paperunits','inches');
set(1,'papersize',[1080 1920]);
set(1,'paperposition',[0 0 1920 1080]);

c_map = gray(4096);

ones_aux_vec = ones(length(ast_ind),1);

for k = source_list 
    
    tic;
    if system_setting_index == 1
    u_data = zeros(floor(length(t_vec)/data_param),n_p);
    f_data = zeros(floor(length(t_vec)/data_param),n_p);
    amp_u = zeros(1,n_p);
    amp_f = zeros(1,n_p);
    end

    rec_data = zeros(floor(length(t_vec)/data_param),n_r);
    t_data = zeros(1,floor(length(t_vec)/data_param));

    data_ind = 0;
    amp_rec = zeros(1,n_r);
    d_ind_aux = d_ind_aux_vec(k);
    if use_gpu
    u = gpuArray(zeros(n_nodes,1));
    p_1 = gpuArray(zeros(n_triangles,1));
    p_2 = gpuArray(zeros(n_triangles,1));
    aux_vec_init_1 = gpuArray(zeros(size(u)));
    aux_vec_init_2 = gpuArray(zeros(size(u)));
    else
    u = zeros(n_nodes,1);
    p_1 = zeros(n_triangles,1);
    p_2 = zeros(n_triangles,1);
    end
    
    
for i = 1 : length(t_vec);

t = (i-1)*d_t;    
f = zeros(n_nodes,1);

if t <= pulse_length
    f(d_ind_aux) = bh_window(t, pulse_length, carrier_freq, 'complex');
end

aux_vec_1 =  - B_1_T*p_1 - B_2_T*p_2 + f;
aux_vec_1 = aux_vec_1 - R*u;
if use_gpu
[aux_vec_2] = pcg_iteration_gpu(C,aux_vec_1,pcg_tol,pcg_max_it,M,aux_vec_init_1); 
aux_vec_init_1 = gpuArray(aux_vec_2);
aux_vec_1 = gpuArray(zeros(size(u)));
else
aux_vec_2 = L_mat\aux_vec_1(p_vec);
aux_vec_2 = U_mat\aux_vec_2;
aux_vec_2 = aux_vec_2(i_p_vec);
aux_vec_1 = zeros(size(u));
end
aux_vec_1(I_p_x) = pml_val*u(I_p_x);
aux_vec_1(I_p_y) = pml_val*u(I_p_y);
u = u + d_t*(aux_vec_2-aux_vec_1);

aux_vec_1 = B_1*u;
if use_gpu 
aux_vec_2 = A.*aux_vec_1;
aux_vec_1 = gpuArray(zeros(size(p_1)));
else
aux_vec_2 = A\aux_vec_1;
aux_vec_1 = zeros(size(p_1));
end
aux_vec_1(I_t_x) = pml_val*p_1(I_t_x);
p_1 = p_1 + d_t*(aux_vec_2 - aux_vec_1);

aux_vec_1 = B_2*u;
if use_gpu 
aux_vec_2 = A.*aux_vec_1;
aux_vec_1 = gpuArray(zeros(size(p_2)));
else
aux_vec_2 = A\aux_vec_1;
aux_vec_1 = zeros(size(p_2));
end
aux_vec_1(I_t_y) = pml_val*p_2(I_t_y);
p_2 = p_2 + d_t*(aux_vec_2 - aux_vec_1);

if mod(i-1,data_param)==0
data_ind = data_ind + 1
aux_vec_2 = B_1_T*p_1 + B_2_T*p_2;
aux_vec_2 = aux_vec_2 + R*u;
if use_gpu 
[aux_vec_2] = pcg_iteration_gpu(C,aux_vec_2,pcg_tol,pcg_max_it,M,aux_vec_init_2); 
aux_vec_init_2 = gpuArray(aux_vec_2);
else
aux_vec_2 = L_mat\aux_vec_2(p_vec);
aux_vec_2 = U_mat\aux_vec_2;
aux_vec_2 = aux_vec_2(i_p_vec);
end
if use_gpu
[rec_data(data_ind,:)] = gather(u(r_ind));
t_data(data_ind) = t;
if system_setting_index == 1
f_data(data_ind,:) = gather(aux_vec_2(ast_p_ind));
u_data(data_ind,:) = gather(u(ast_p_ind));
end
else
rec_data(data_ind,:) = u(r_ind);
t_data(data_ind) = t;
if system_setting_index == 1    
f_data(data_ind,:) = aux_vec_2(ast_p_ind);
u_data(data_ind,:) = u(ast_p_ind);
end
end
end

if i == frame_param
figure(1); clf; hold on; set(gcf,'renderer','opengl'); set(gcf,'color',[1 1 1]);
u_plot = gather(real(u)/max(abs(real(u(:)))));
h_surf=trisurf(triangles(:,1:3),nodes(:,1),nodes(:,2),zeros(size(nodes,1),1),u_plot);
view(0,90); set(gca,'visible','off');
shading interp; lighting phong; camlight right; camlight left; material dull; 
plot_boundary_curve;

colormap(c_map);
set(gca,'clim',[-1 1]);
set(gca,'xlim',[-0.27 0.27]);
set(gca,'ylim',[-0.27 0.27]);
axis equal;
drawnow;

elseif mod(i,frame_param)==0
u_plot = gather((real(u))/max(abs(real(u(:)))));
set(h_surf,'cdata',u_plot);
drawnow;

end

end

if system_setting_index == 1
save(['./' data_name '/f_data_modulated_' int2str(k) '.mat'],'f_data','t_data');
save(['./' data_name '/u_data_modulated_' int2str(k) '.mat'],'u_data','rec_data','t_data');
else
save(['./' data_name '/u_data_modulated_' int2str(k) '.mat'],'rec_data');
end

toc

end
