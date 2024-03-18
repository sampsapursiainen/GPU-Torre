%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

load data/triangles_2.dat
load data/nodes_2.dat

system_setting_index = 2;
parameters;

triangles = triangles_2;
nodes = nodes_2;

triangles = [triangles zeros(size(triangles,1),3)];

edges_aux = [triangles(:,[2 3]) [1:size(triangles,1)]'  1*ones(size(triangles,1),1);
             triangles(:,[3 1]) [1:size(triangles,1)]'  2*ones(size(triangles,1),1);
             triangles(:,[1 2]) [1:size(triangles,1)]'  3*ones(size(triangles,1),1)];
[edges, edge_ind_1, edge_ind_2] = unique(sort(edges_aux(:,1:2),2),'rows');

triangles(4*size(triangles,1)+(edges_aux(:,4)-1)*size(triangles,1)+edges_aux(:,3)) = edge_ind_2;

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
cover_ind = find(ismember(triangles(:,4),ast_aux_ind(2)));
ast_ind_1 = find(ismember(triangles(:,4),ast_aux_ind(1)));

N_s = 150;
z_s = rand(N_s,N_s);
z_s = z_s - mean(z_s(:));
W_s = speye(N_s^2,N_s^2);
phi_s = pi/4;

Ind_mat_s = reshape([1 : N_s^2]',N_s,N_s); 

W_s = 10*speye(N_s^2,N_s^2);

for i = 2 : N_s-1
    for j = 2 : N_s-1 
            W_s(Ind_mat_s(i,j),Ind_mat_s(i,j+1)) = W_s(Ind_mat_s(i,j),Ind_mat_s(i,j+1)) - 1;
            W_s(Ind_mat_s(i,j),Ind_mat_s(i-1,j)) = W_s(Ind_mat_s(i,j),Ind_mat_s(i-1,j)) - 1;
            W_s(Ind_mat_s(i,j),Ind_mat_s(i,j-1)) = W_s(Ind_mat_s(i,j),Ind_mat_s(i,j-1)) - 1;
            W_s(Ind_mat_s(i,j),Ind_mat_s(i+1,j)) = W_s(Ind_mat_s(i,j),Ind_mat_s(i+1,j)) - 1;
    end
end
x_s = W_s \ z_s(:);
x_s = x_s - min(x_s) + 1e-6; 
x_s = x_s / max(abs(x_s));

nodes_c = (nodes(triangles(ast_ind,1),:) +  nodes(triangles(ast_ind,2),:) +  nodes(triangles(ast_ind,3),:))/3;
I_1 = floor((N_s-1)*(nodes_c(:,1)+0.16)/0.32 + 1);
I_2 = floor((N_s-1)*(-nodes_c(:,2)+0.16)/0.32 + 1);
x_s = x_s(N_s*(I_1-1) + I_2);
I_1 = find(x_s>0.5);
I_2 = find(x_s<0.5);
x_s(I_1) = 1;
x_s(I_2) = 0;

p_vec(ast_ind) = bg_permittivity;
p_vec(cover_ind) = sl_permittivity;

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
D = spalloc(n_nodes,n_nodes,0);

ind_e_1 = [2 3;
           3 1;
           1 2];
ind_e_2 = [ 0  3  -2; 
           -3  0  1;
           2 -1  0];
       
t_c_vec = (1/3)*(nodes(triangles(:,1),:) + nodes(triangles(:,2),:) + nodes(triangles(:,3),:));
I_t_x = find(abs(t_c_vec(:,1))>0.2);
I_t_y = find(abs(t_c_vec(:,2))>0.2);
I_p_x = unique(triangles(I_t_x,1:3));
I_p_y = unique(triangles(I_t_y,1:3));
I_ast = unique(triangles(ast_ind,1:3));
I_or = find(sqrt(sum(t_c_vec.^2,2)) > or_radius);
I_or= unique(triangles(I_or,1:3));

save data/mesh_2.mat n_nodes n_triangles n_edges triangles nodes edges ind_e_1 ind_e_2 p_vec I_t_x I_t_y I_p_x I_p_y I_ast I_or interp_vec;       
       
Aux_mat = [nodes(triangles(:,1),:)' ; nodes(triangles(:,2),:)'] - repmat(nodes(triangles(:,3),:)',2,1); 
ind_m = [1 3; 2 4];
ala = abs((Aux_mat(ind_m(1,1),:).*Aux_mat(ind_m(2,2),:)-Aux_mat(ind_m(1,2),:).*Aux_mat(ind_m(2,1),:)))/2;
clear Aux_mat; 

grad_arr = zeros(6,size(triangles,1));

for i = 1 : 3
 
grad_arr(2*(i-1)+1:2*i,:) = nodes(triangles(:,ind_e_1(i,2)),[2 1])'-nodes(triangles(:,ind_e_1(i,1)),[2 1])';  
grad_arr(2*(i-1)+1,:) = - grad_arr(2*(i-1)+1,:);
aux_vec =  dot(grad_arr(2*(i-1)+1:2*i,:),(nodes(triangles(:,i),:)'-nodes(triangles(:,ind_e_1(i,1)),:)'));
grad_arr(2*(i-1)+1:2*i,:) = (1./aux_vec([1 1],:)).*grad_arr(2*(i-1)+1:2*i,:);   

end

for i = 1 : 3
    %i
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
    %i
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


for i = 1 : 3
    %i
    for j = i : 3
        
         d_vec = zeros(1,size(triangles,1));
                 
           if i == j
            aux_val = 1/6;
           else
            aux_val = 1/12;
           end
          
           d_vec = d_vec + aux_val.*ala;
                          
       D_part = sparse(triangles(:,i), triangles(:,j), d_vec',n_nodes,n_nodes);
        
       if i == j
           D = D + D_part;
       else
           D = D + D_part;
           D = D + D_part';
       end
       
    end
end


A = spdiags(ala',0,n_triangles,n_triangles);


for i = 1 : 3
    %i   
    
       aux_vec = grad_arr(2*(i-1)+1:2*i,:);
       b_vec = ala([1 1],:).*aux_vec;
                          
       B_1_part = sparse([1:size(triangles,1)]', triangles(:,i), b_vec(1,:)',n_triangles,n_nodes);
       B_2_part = sparse([1:size(triangles,1)]', triangles(:,i), b_vec(2,:)',n_triangles,n_nodes); 
       
       B_1 = B_1 + B_1_part;
       B_2 = B_2 + B_2_part;       
   
end


int_vec = zeros(n_nodes,1);

for i = 1 : 3
    %i
                 
       aux_val = 1/3; 
         
       aux_vec = aux_val.*ala';
                          
       int_vec(triangles(:,i)) = int_vec(triangles(:,i)) + aux_vec;
        

end

J_aux_1 = sparse(size(C,1),length(ast_ind));
J_aux_2 = sparse(size(C,1),length(ast_ind));
J_aux_3 = sparse(size(C,1),length(ast_ind));



save data/system_data_2.mat A B_1 B_2 R C D J_aux_1 J_aux_2 J_aux_3 r_ind rec_pos ast_ind int_vec;