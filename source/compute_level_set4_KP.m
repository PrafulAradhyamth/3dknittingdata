function [level_set,tri_idx,level_field] = compute_level_set4_KP(TR,edges_unique, level,D,scalar_field)
% compute all points on triangulation edges with a specified value of the
% given scalar field, i.e. compute points on a isocontour
% input:
% TR: triangulation
% edges_unique: triangulation edges
% level: specified value of the scalar field
% D: scalar field in question
% scalar field: cell of possible other scalar fields, i.e. distance from
% middle line
% 
% output:
% level_set: coordinates in the level set; they lay on triangle edges
% tri_idx: inidces of triangles which contain edges with level_set points
% (important for sorting)
% level_field: values of other scalar field at level_set
    
neg_tri = -1;
level_set = [];
level_set_points = [];
tri_idx = [];
level_field = cell(1,length(scalar_field));
% just in case, snap distances with small abs value to 0
D(abs(D)<1e-3) = 0;

D_edges = D(edges_unique);
vertex = TR.Points';
% check all edges. 
% if one vertex values is below and the other one above the level,
% there is a contour point somewhere on the edge
for i=1:length(edges_unique)
    if prod(D_edges(i,:)-level) <= 0 % if one edge vertex is above and the other below the level
        [D_sorted,idx] = sort(D_edges(i,:));
        vtx1 = edges_unique(i,idx(1));
        vtx2 = edges_unique(i,idx(2));
        
        p1 = vertex(:,vtx1);
        p2 = vertex(:,vtx2);
        t0 = (level-D_sorted(1))/(D_sorted(2)-D_sorted(1));
        if isnan(t0); continue; end
        p0 = (1-t0)*p1 + t0*p2; % compute the point on the edge 
        if p0 == p1 | p0 == p2
            level_set_points = cat(2,level_set_points,p0);
        end
        for j=1:length(scalar_field) % compute the corresponing value for possible
            % other scalar fields
            field1 = scalar_field{j}(vtx1);
            field2 = scalar_field{j}(vtx2);
            field0 = (1-t0)*field1 + t0*field2;
            level_field{j} = cat(1,level_field{j},field0);
        end

        level_set = cat(2,level_set,p0);
        % which triangles belong to edges with points on the contour 
        % (important for sorting the contour points)
        r = edgeAttachments(TR,edges_unique(i,idx(1)),edges_unique(i,idx(2)));
        % boundary 
        r = sort(r{1})';
        if length(r) < 2
            r = [r(1);neg_tri];
            neg_tri = neg_tri-1;
        end
        tri_idx = cat(1,tri_idx,r');
    end
end
if not(isempty(level_set_points))
    %Bereinigung wenn Punkte in den Ecken der Kanten liegen, bzw. Punkte
    %mehrere relevante Kanten hervorruven
   [tri_idx, level_set, level_field] = findContourPathKP(tri_idx,level_set, level_field, TR);
end
[tri_idx,sidx] = sortrows(tri_idx);
level_set = level_set(:,sidx);
for j=1:length(level_field)
    level_field{j} = level_field{j}(sidx);
end