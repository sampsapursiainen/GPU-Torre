%Copyright © 2019- Sampsa Pursiainen & GPU-ToRRe Development Team
%See: https://github.com/sampsapursiainen/GPU-Torre

g = decsg(gd,sf,ns);

h_max = 0.08;
[p,e,t]=initmesh(g,'hmax',h_max,'hgrad',1.08);
p = jigglemesh(p,e,t,'iter',1);

figure(1); clf; pdeplot(p,e,t);
set(gca,'xlim',[-0.75 0.75]);
set(gca,'ylim',[-0.75 0.75]);
axis equal;
n_triangles = size(t,2)
n_nodes = size(p,2)

t = t';
p = p';

save data/triangles_1.dat -ascii t;
save data/nodes_1.dat -ascii p;
