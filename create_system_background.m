%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
load data/triangles_1.dat
load data/nodes_1.dat

i_ind = 1;
parameters;

triangles = triangles_1;
nodes = nodes_1;

triangles = [triangles zeros(size(triangles,1),3) ];

edges_aux = [triangles(:,[2 3]) [1:size(triangles,1)]'  1*ones(size(triangles,1),1);
             triangles(:,[3 1]) [1:size(triangles,1)]'  2*ones(size(triangles,1),1);
             triangles(:,[1 2]) [1:size(triangles,1)]'  3*ones(size(triangles,1),1)];
[edges, edge_ind_1, edge_ind_2] = unique(sort(edges_aux(:,1:2),2),'rows');

triangles(4*size(triangles,1)+(edges_aux(:,4)-1)*size(triangles,1)+edges_aux(:,3)) = edge_ind_2;

ast_ind_aux = find(triangles(:,4)==ast_aux_ind);
nodes_aux = (nodes(triangles(ast_ind_aux,1),:) + nodes(triangles(ast_ind_aux,2),:) + nodes(triangles(ast_ind_aux,3),:))/3;
n_nodes_aux = size(nodes_aux,1);
tri_aux = delaunay(nodes_aux(:,1), nodes_aux(:,2));
edges_aux = [tri_aux(:,[2 3]) [1:size(tri_aux,1)]'  1*ones(size(tri_aux,1),1);
             tri_aux(:,[3 1]) [1:size(tri_aux,1)]'  2*ones(size(tri_aux,1),1);
             tri_aux(:,[1 2]) [1:size(tri_aux,1)]'  3*ones(size(tri_aux,1),1)];
edges_aux = sortrows(sort(edges_aux(:,1:2),2));

I_boundary_aux = find(sum(abs(edges_aux - [0 0 ; edges_aux(1:end-1,:)]),2)>0 & sum(abs(edges_aux - [edges_aux(2:end,:); 0 0]),2)>0);
I_boundary_aux = unique(edges_aux(I_boundary_aux,:));

G_aux = spalloc(n_nodes_aux, n_nodes_aux,0);

Aux_mat = [nodes_aux(tri_aux(:,1),:)' ; nodes_aux(tri_aux(:,2),:)'] - repmat(nodes_aux(tri_aux(:,3),:)',2,1); 
ind_m = [1 3; 2 4];
ala_aux = abs((Aux_mat(ind_m(1,1),:).*Aux_mat(ind_m(2,2),:)-Aux_mat(ind_m(1,2),:).*Aux_mat(ind_m(2,1),:)))/2;

ind_e_1 = [2 3;
           3 1;
           1 2];
ind_e_2 = [ 0  3  -2; 
           -3  0  1;
           2 -1  0];

grad_arr_aux = zeros(6,size(tri_aux,1));

for i = 1 : 3
 
grad_arr_aux(2*(i-1)+1:2*i,:) = nodes_aux(tri_aux(:,ind_e_1(i,2)),[2 1])'-nodes_aux(tri_aux(:,ind_e_1(i,1)),[2 1])';  
grad_arr_aux(2*(i-1)+1,:) = - grad_arr_aux(2*(i-1)+1,:);
aux_vec = dot(grad_arr_aux(2*(i-1)+1:2*i,:),(nodes_aux(tri_aux(:,i),:)'-nodes_aux(tri_aux(:,ind_e_1(i,1)),:)'));
grad_arr_aux(2*(i-1)+1:2*i,:) = (1./aux_vec([1 1],:)).*grad_arr_aux(2*(i-1)+1:2*i,:);   

end

for i = 1 : 3
   % i
    for j = i : 3
        
         g_vec = zeros(1,size(tri_aux,1));
                 
          
           g_vec = g_vec + dot(grad_arr_aux(2*(i-1)+1:2*i,:),grad_arr_aux(2*(j-1)+1:2*j,:)).*ala_aux;
                          
       G_part = sparse(tri_aux(:,i), tri_aux(:,j), g_vec',n_nodes_aux,n_nodes_aux);
        
       if i == j
           G_aux = G_aux + G_part;
       else
           G_aux = G_aux + G_part;
           G_aux = G_aux + G_part';
       end
       
    end
end

    ast_p_ind = unique(triangles(ast_ind_aux,1:3));
    triangles_ast = triangles(ast_ind_aux,1:3);
    aux_p_ind = zeros(size(nodes,1),1);
    aux_p_ind(ast_p_ind) = [1:length(ast_p_ind)]';
    triangles_ast = aux_p_ind(triangles_ast);
    
interp_vec = [1:size(triangles,1)];
for n_ref_aux_ind = 1 : n_ref_aux
[nodes,edges,triangles,interp_vec_aux] = mesh_refinement(nodes,edges,triangles);
interp_vec = interp_vec(interp_vec_aux);
end

r_phi = [0:2*pi/n_r:2*pi*(1-1/n_r)]';
rec_pos = or_radius*[cos(r_phi) sin(r_phi)];
r_ind = zeros(1,n_r);
for i = 1 : n_r
[r_val_aux, r_ind_val] = min(sum((nodes - repmat(rec_pos(i,:),size(nodes,1),1)).^2,2));
r_ind(i) = r_ind_val;
end

p_vec = ones(size(triangles,1),1);
ast_ind = find(ismember(triangles(:,4),ast_aux_ind));



p_vec(ast_ind) = bg_permittivity;
for i = 1 : length(cav_ind)
I = find(ismember(triangles(:,4),cav_ind(i)));
p_vec(I) = cav_val(i);
end

n_nodes = size(nodes,1);
n_triangles = size(triangles,1);
n_edges = size(edges,1);

C = spalloc(n_nodes,n_nodes,0);
R = spalloc(n_nodes,n_nodes,0);
A = spalloc(n_triangles,n_triangles,0);
B_1 = spalloc(n_triangles,n_nodes,0);
B_2 = spalloc(n_triangles,n_nodes,0);
D = spalloc(length(ast_p_ind),length(ast_p_ind),0);
       
t_c_vec = (1/3)*(nodes(triangles(:,1),:) + nodes(triangles(:,2),:) + nodes(triangles(:,3),:));
I_t_x = find(abs(t_c_vec(:,1))>0.2);
I_t_y = find(abs(t_c_vec(:,2))>0.2);
I_p_x = unique(triangles(I_t_x,1:3));
I_p_y = unique(triangles(I_t_y,1:3));
I_ast = unique(triangles(ast_ind,1:3));
I_or = find(sqrt(sum(t_c_vec.^2,2)) > or_radius);
I_or= unique(triangles(I_or,1:3));
%I_b_y = setdiff(I_b_y,I_b_x);

save data/mesh_1.mat n_nodes n_triangles n_edges triangles nodes edges ind_e_1 ind_e_2 p_vec I_t_x I_t_y I_p_x I_p_y I_ast I_or interp_vec;       
       
Aux_mat = [nodes(triangles(:,1),:)' ; nodes(triangles(:,2),:)'] - repmat(nodes(triangles(:,3),:)',2,1); 
Aux_mat_1 = [nodes(triangles(ast_ind,1),:)' ; nodes(triangles(ast_ind,2),:)'] - repmat(nodes(triangles(ast_ind,3),:)',2,1); 
ind_m = [1 3; 2 4];
ala = abs((Aux_mat(ind_m(1,1),:).*Aux_mat(ind_m(2,2),:)-Aux_mat(ind_m(1,2),:).*Aux_mat(ind_m(2,1),:)))/2;
ala_1 = abs((Aux_mat_1(ind_m(1,1),:).*Aux_mat_1(ind_m(2,2),:)-Aux_mat_1(ind_m(1,2),:).*Aux_mat_1(ind_m(2,1),:)))/2;
clear Aux_mat; 

grad_arr = zeros(6,size(triangles,1));

for i = 1 : 3
 
grad_arr(2*(i-1)+1:2*i,:) = nodes(triangles(:,ind_e_1(i,2)),[2 1])'-nodes(triangles(:,ind_e_1(i,1)),[2 1])';  
grad_arr(2*(i-1)+1,:) = - grad_arr(2*(i-1)+1,:);
aux_vec =  dot(grad_arr(2*(i-1)+1:2*i,:),(nodes(triangles(:,i),:)'-nodes(triangles(:,ind_e_1(i,1)),:)'));
grad_arr(2*(i-1)+1:2*i,:) = (1./aux_vec([1 1],:)).*grad_arr(2*(i-1)+1:2*i,:);   

end


for i = 1 : 3
  %  i
    for j = i : 3
        
         c_vec = zeros(1,size(triangles,1));
                 
           if i == j
            aux_val = 1/6;
           else
            aux_val = 1/12;
           end
          
       c_vec = c_vec + absorb_val*p_vec'.*aux_val.*ala;
                          
       C_part = sparse(triangles(:,i), triangles(:,j), c_vec',n_nodes,n_nodes);
        
       if i == j
           R = R + C_part;
       else
           R = R + C_part;
           R = R + C_part';
       end
       
    end
end

for i = 1 : 3
   % i
    for j = i : 3
        
         c_vec = zeros(1,size(triangles,1));
                 
           if i == j
            aux_val = 1/6;
           else
            aux_val = 1/12;
           end
          
           c_vec = c_vec + p_vec'.*aux_val.*ala;
                          
       C_part = sparse(triangles(:,i), triangles(:,j), c_vec',n_nodes,n_nodes);
        
       if i == j
           C = C + C_part;
       else
           C = C + C_part;
           C = C + C_part';
       end
       
    end
end

A = spdiags(ala',0,n_triangles,n_triangles);

for i = 1 : 3
 %   i   
    
       aux_vec = grad_arr(2*(i-1)+1:2*i,:);
       b_vec = ala([1 1],:).*aux_vec;
                          
       B_1_part = sparse([1:size(triangles,1)]', triangles(:,i), b_vec(1,:)',n_triangles,n_nodes);
       B_2_part = sparse([1:size(triangles,1)]', triangles(:,i), b_vec(2,:)',n_triangles,n_nodes); 
       
       B_1 = B_1 + B_1_part;
       B_2 = B_2 + B_2_part;       
   
end


int_vec = zeros(n_nodes,1);

for i = 1 : 3
%    i
                 
       aux_val = 1/3; 
         
       aux_vec = aux_val.*ala';
                          
       int_vec(triangles(:,i)) = int_vec(triangles(:,i)) + aux_vec;
        
end

save data/system_data_1.mat A B_1 B_2 R C G_aux r_ind rec_pos ast_ind ast_p_ind triangles_ast int_vec nodes_aux tri_aux;