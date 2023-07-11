function [seam,TR_new,corresponding_vertices_ID] = insert_seam(sm,TR,edges_unique)
    % insert seam and middle line into the mesh:
    % input: 
    % sm: matrix of points
    % TR: triangulation
    % edges_unique: triangulation edges
    % output:
    % seam: (x,y,z) coordinates of seam/mid points
    % TR_new: updated triangulation with inserted points
    % corresponding_vetrtices_ID: indices of new points 
    
    
    % sm has the following form: x,y,z,index,type
    % type can be:
    %   - vertex: 1
    %   - edge:   2  
    %   - face:   3
    seam = sm(:,1:3);
    conl = TR.ConnectivityList;
    vtx = TR.Points;
    corresponding_vertices_ID = zeros(size(sm,1),1);
    for i=1:size(sm,1)
        if sm(i,5) == 1
            % if an existing point, just note
            corresponding_vertices_ID(i) = sm(i,4);
        elseif sm(i,5) == 3
            error('Inserting vertex on faces is not supported!')
        else
            % if on edge, insert point into the mesh
            vtx_idx = edges_unique(sm(i,4),:);
            
            
            [conl,vtx,new_vertex_ID,~] = insert_vertex(conl,vtx,seam(i,:),vtx_idx);
            corresponding_vertices_ID(i) = new_vertex_ID;
        end       
    end
    TR_new = triangulation(conl,vtx);