%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
function [nodes,edges,triangles,interp_vec] = mesh_refinement(nodes,edges,triangles)

n_nodes = size(nodes,1);
n_triangles = size(triangles,1);

nodes = [nodes;
         (1/2)*(nodes(edges(:,1),:) + nodes(edges(:,2),:))];

interp_vec = [1:n_triangles]';    
interp_vec = interp_vec(:,[1 1 1 1]);
interp_vec = interp_vec(:);

t_aux_1 = triangles(:,1);
t_aux_2 = triangles(:,2);
t_aux_3 = triangles(:,3);
t_aux_4 = triangles(:,4);
t_aux_5 = n_nodes+triangles(:,5);
t_aux_6 = n_nodes+triangles(:,6);
t_aux_7 = n_nodes+triangles(:,7);

triangles = [t_aux_1 t_aux_7 t_aux_6 t_aux_4; 
             t_aux_7 t_aux_2 t_aux_5 t_aux_4;
             t_aux_5 t_aux_3 t_aux_6 t_aux_4;
             t_aux_5 t_aux_6 t_aux_7 t_aux_4];

triangles = [triangles zeros(size(triangles,1),3)];

edges_aux = [triangles(:,[2 3]) [1:size(triangles,1)]'  1*ones(size(triangles,1),1);
             triangles(:,[3 1]) [1:size(triangles,1)]'  2*ones(size(triangles,1),1);
             triangles(:,[1 2]) [1:size(triangles,1)]'  3*ones(size(triangles,1),1)];
[edges, edge_ind_1, edge_ind_2] = unique(sort(edges_aux(:,1:2),2),'rows');

triangles(4*size(triangles,1)+(edges_aux(:,4)-1)*size(triangles,1)+edges_aux(:,3)) = edge_ind_2;