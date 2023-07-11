function [pt_order,eu] = contour_forming2(nodes,H,valb)
    % ordering of triangle indices
    % input:
    % nodes: connected components of graph representing the ordering of
    % triangle edges
    % H: the whole graph
    % valb: indices of triangle edges in H
    % output:
    % pt_order: ordered point (i.e. triangle) indices
    % eu: label if an Eulerian path exists (see grIsEulerian)
    
Hc = subgraph(H,nodes); 
EdgesHc = nodes(table2array(Hc.Edges));
[eu,cEu] = grIsEulerian(table2array(Hc.Edges)); 
concomEdges = EdgesHc(cEu,:);
val_Hc = zeros(size(concomEdges,1),1);
for k=1:size(concomEdges,1) 
    val_Hc(k) = findedge(H,concomEdges(k,1),concomEdges(k,2)); %reordered edges according to cEu
end
reorder_Hc = zeros(length(val_Hc),1);
clear k
for k=1:length(val_Hc); reorder_Hc(k) = find(valb == val_Hc(k)); end
if eu == 1
    reorder_Hc = [reorder_Hc(:); reorder_Hc(1)];
end
pt_order = reorder_Hc;