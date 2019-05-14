%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
function [path_data] = create_constellation(n_r, m, bistatic_difference, orbiting_velocity, sparsity_factor, back_scattering,plot_constellation,case_num) 

ind_aux = [1:sparsity_factor:n_r]';
path_data = ind_aux;
path_data_aux = [];
aux_ind = [0:bistatic_difference*sparsity_factor:n_r];
aux_ind = aux_ind(2:m);
for i = 1 : m-1
aux_vec = mod(floor(aux_ind(i) + orbiting_velocity*ind_aux-1), n_r) + 1;
path_data = [path_data aux_vec];
path_data_aux = [path_data_aux ; ind_aux aux_vec];  
end
if back_scattering == 1
aux_vec = ind_aux;
path_data = [path_data aux_vec];
path_data_aux = [path_data_aux ; ind_aux aux_vec];  
end

if plot_constellation
X = [cos([0 : 2*pi/n_r : (n_r-1)*2*pi/n_r]); sin([0 : 2*pi/n_r : (n_r - 1)*2*pi/n_r])];
fig = figure(2); 
clf; 
h_plot = plot([X(1,path_data_aux(:,1)) ;X(1,path_data_aux(:,2))],[X(2,path_data_aux(:,1)) ;X(2,path_data_aux(:,2))]); 
set(h_plot,'linewidth',1);
set(h_plot,'color','k');
axis equal; 
set(gca,'visible','off'); 
set(gcf,'renderer','zbuffer');
hold on;
t_vec = [0 : 0.01 : 2*pi];
plot(cos(t_vec),sin(t_vec),'k','linewidth',1)
set(gca,'xlim',[-1.01 1.01])
set(gca,'ylim',[-1.01 1.01])
drawnow;
hold off;
end

print(fig,['constellation_' case_num],'-dpng')

path_data_aux = [];
for i = 1 : size(path_data,2)-1
    path_data_aux = [path_data_aux ; [path_data(:,1) path_data(:,i+1)]];
end
path_data = path_data_aux;

end
