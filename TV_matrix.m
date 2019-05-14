%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre
function TV_D = TV_matrix(nodes,triangles)

edges_aux = [triangles(:,[2 3]) [1:size(triangles,1)]' 
             triangles(:,[3 1]) [1:size(triangles,1)]'  
             triangles(:,[1 2]) [1:size(triangles,1)]' ];
edges_aux(:,1:2) = sort(edges_aux(:,1:2),2);
edges_aux = sortrows(edges_aux);
I = find(sum((edges_aux(1:end-1,1:2)-edges_aux(2:end,1:2)).^2,2)==0);
edges_aux = [edges_aux(I,:) ; edges_aux(I+1,:)]; 
edges_aux = sortrows(edges_aux);

TV_D = sparse(size(triangles,1), size(triangles,1)); 
D_aux = sparse(size(nodes,1), size(nodes,1));

edge_length = sqrt(sum((nodes(edges_aux(:,1),:) - nodes(edges_aux(:,2),:)).^2,2));

for i = 1 : length(edge_length)

TV_D(edges_aux(i,3), edges_aux(i,3)) = TV_D(edges_aux(i,3), edges_aux(i,3)) + edge_length(i);    

if D_aux(edges_aux(i,1), edges_aux(i,2)) == 0
D_aux(edges_aux(i,1), edges_aux(i,2)) = edges_aux(i,3);
D_aux(edges_aux(i,2), edges_aux(i,1)) = edges_aux(i,3);
elseif D_aux(edges_aux(i,1), edges_aux(i,2)) > 0 
 aux_val = D_aux(edges_aux(i,1), edges_aux(i,2));
 TV_D(edges_aux(i,3), aux_val) = TV_D(edges_aux(i,3), aux_val) - edge_length(i);
 TV_D(aux_val, edges_aux(i,3)) = TV_D(aux_val, edges_aux(i,3)) - edge_length(i);
end

end

TV_D = full(TV_D)/max(edge_length);

end