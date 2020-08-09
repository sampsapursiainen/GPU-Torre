%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

load data/mesh_1.mat;
load data/system_data_1.mat; 

fig = figure(2); clf;
c_map = gray(4096);
c_map = flipud(c_map);
plot_x_coords = nodes(ast_p_ind(triangles_ast'),1);
plot_y_coords = nodes(ast_p_ind(triangles_ast'),2);
plot_x_coords = reshape(plot_x_coords,3,length(plot_x_coords)/3);
plot_y_coords = reshape(plot_y_coords,3,length(plot_y_coords)/3);
h_ast = patch(plot_x_coords, plot_y_coords, max(db(abs(x)/max(abs(x))),plot_threshold_db));
set(h_ast,'edgecolor','none');
colormap(c_map);
set(gca,'visible','off');
colorbar;
axis equal;
%set(gca,'clim',[-3 1]);
