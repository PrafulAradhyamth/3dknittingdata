function TR = remedy_mesh(TR)
% remove all copies of multiple triangles
F = TR.ConnectivityList;
V = TR.Points;
TRF_sorted = sort(F,2);
[~,~,ic] = unique(TRF_sorted, 'rows', 'stable');        
h = accumarray(ic, 1);                             
maph = h(ic);
tr_to_delete = maph>1;

F(tr_to_delete,:) = [];

[V2, F2] = trimMesh(V, F);
TR = triangulation(F2,V2);