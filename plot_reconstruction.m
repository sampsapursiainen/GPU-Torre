%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
load data/mesh_1.mat;
load data/system_data_1.mat; 

fig = figure(2); clf;
plot_x_coords = nodes(ast_p_ind(triangles_ast'),1);
plot_y_coords = nodes(ast_p_ind(triangles_ast'),2);
plot_x_coords = reshape(plot_x_coords,3,length(plot_x_coords)/3);
plot_y_coords = reshape(plot_y_coords,3,length(plot_y_coords)/3);
h_ast = patch(plot_x_coords, plot_y_coords, x + 4);
set(h_ast,'edgecolor','none');
colormap gray;
set(gca,'visible','off');
colorbar;
axis equal;
set(gca,'clim',[1 5]);
