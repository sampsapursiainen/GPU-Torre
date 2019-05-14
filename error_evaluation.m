%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
plot_reconstruction;

domain_size = 0.15;
im_res = 1024;
[X_im,Y_im] = meshgrid(linspace(-domain_size,domain_size,im_res));

load data/mesh_2.mat;
load data/system_data_2.mat;

c_tri = (1/3)*(nodes(triangles(:,1),:) + nodes(triangles(:,2),:) + nodes(triangles(:,3),:));

mask_surface = zeros(size(triangles,1),1);
mask_surface(find(triangles(:,4) == 3)) = 1;

mask_void = zeros(size(triangles,1),1);
mask_void(find(ismember(triangles(:,4),[5,6,7]))) = 1;

mask_asteroid = zeros(size(triangles,1),1);
mask_asteroid(find(ismember(triangles(:,4),[2,3,5,6,7]))) = 1;

mask_void = griddata(c_tri(:,1),c_tri(:,2),mask_void,X_im,Y_im,'nearest');
mask_surface = griddata(c_tri(:,1),c_tri(:,2),mask_surface,X_im,Y_im,'nearest');
mask_asteroid = griddata(c_tri(:,1),c_tri(:,2),mask_asteroid,X_im,Y_im,'nearest');

P_im_exact = griddata(c_tri(:,1),c_tri(:,2),p_vec,X_im,Y_im,'nearest');

%figure(3);
%clf;
%pcolor(P_im_exact);
%shading flat;
%set(gca,'visible','off');
%colorbar;
%colormap gray;

load data/triangles_1.dat;
load data/nodes_1.dat;

ast_ind = find(triangles_1(:,4)==3);

c_tri = (1/3)*(nodes_1(triangles_1(:,1),:) + nodes_1(triangles_1(:,2),:) + nodes_1(triangles_1(:,3),:));

p_vec = ones(size(triangles_1,1),1);
p_vec(ast_ind) = x + 4;

P_im_reconstruction = griddata(c_tri(:,1),c_tri(:,2),p_vec,X_im,Y_im,'nearest');

%figure(4);
%clf;
%pcolor(P_im_reconstruction);
%shading flat;
%set(gca,'visible','off');
%colorbar;
%colormap gray;

error_global = immse(P_im_exact,P_im_reconstruction)
error_void = immse(mask_void.*P_im_exact,mask_void.*P_im_reconstruction)
error_surface = immse(mask_surface.*P_im_exact,mask_surface.*P_im_reconstruction)

aux_length = length(find(mask_void + mask_surface));

aux_ind_1 = find(mask_asteroid);
aux_vec = sort(P_im_reconstruction(aux_ind_1));
thresh_val = aux_vec(aux_length);

P_im_overlap = zeros(size(P_im_reconstruction));
P_im_overlap(find(P_im_reconstruction.*mask_asteroid < thresh_val)) = 1; 

aux_val_1 = length(find(mask_surface));
aux_val_2 = length(find(P_im_overlap.*mask_surface));

figure(3);
clf;
pcolor(-(P_im_overlap+0.5).*mask_surface-mask_asteroid);
shading flat;
set(gca,'visible','off');
colormap gray;
axis equal;

overlap_error_surface = 100 - 100*aux_val_2/aux_val_1

aux_val_1 = length(find(mask_void)); 
aux_val_2 = length(find(P_im_overlap.*mask_void));

figure(4);
clf; 
pcolor(-(P_im_overlap+0.5).*mask_void-mask_asteroid);
shading flat;
set(gca,'visible','off');
colormap gray;
axis equal;

overlap_error_void = 100 - 100*aux_val_2/aux_val_1

