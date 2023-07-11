function ET = edge_triangle_adjacency_DS(F,n,E)
  % EDGE_TRIANGLE_ADJACENCY Build an edge-triangle adjacency matrix for a
  % triangle mesh.
  %
  % ET = edge_triangle_adjacency(F,E)
  %
  % Input:
  %   F  #F by 3  matrix of indices of vertices at triangle corners
  %   E  #E by 2  matrix of indices of vertices at edge endpoints
  % Output:
  %   ET #E by 2  map between an edge to its two indicent faces
  %               (-1 in column 2 if the edge is on the border)
  %
  % Example:
  %  E = edges(F);
  %  ET = edge_triangle_adjacency(F,E);
  % works only for one edge!!!
%   VT = vertex_triangle_adjacency(F);
  [VF,NI] = vertex_triangle_adjacency_mex(F,n);
%   vtxid = find(verticeL);
%   indexfcn = @(s,e)VF(s:e)';
%   faceIDs = arrayfun(indexfcn,NI(vtxid)+1,NI(vtxid+1),'UniformOutput',0)';
%   VT2 = [faceIDs{:}];
  CFi1 = VF(NI(E(1))+1:NI(E(1)+1));
  CFi2 = VF(NI(E(2))+1:NI(E(2)+1));
%   ET = intersect(CFi1,CFi2);
  ET = MY_intersect(CFi1,CFi2);



end
