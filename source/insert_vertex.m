function [faces,vertices,new_vertex_ID,edge_faces] = insert_vertex(faces,vertices,new_vertex,vtx_idx)

    % insert vertex into the mesh
    % input:
    % faces: faces in the mesh (triangles)
    % vertices: vertices in the mesh
    % new_vertex: new vertex coordinates
    % vtx_idx: indices od vertices belonging to the edge where the new vertex lays
    % output:
    % faces: updated faces
    % vertices: updated vertices
    % new_vertex_ID: index of the new vertex
    % edge_faces: index of faces adjacent to all edges


%  cases:
%   - Point is on a boundary edge
%   - Point is on an edge inside the mesh


    edge_verts = vtx_idx;
    
    edge_faces = edge_triangle_adjacency_DS(faces,size(vertices,1),edge_verts);

%     
%     % Two attached faces on edge
        
        % Opposite vertices
        edge_faces_def = faces(edge_faces,:);      % (2 x 3)
        Lia = ismember(edge_faces_def,edge_verts); % logical (2 x 3), true where vertice IDs are on edge 
        replace_face1 = find(Lia(1,:));            % (1 x 2) position of edge's vertice IDs in edge face 1  
        if numel(edge_faces)==2
            replace_face2 = find(Lia(2,:));            % (1 x 2) position of edge's vertice IDs in edge face 2  
        end
        % Add 1 new vertex
        vertices = [vertices;new_vertex];  % (#vertices x 3)
        new_vertex_ID = size(vertices,1);  % scalar integer
% 
        
        % Add 4 new faces
        if numel(edge_faces)==2
            new_faces = [            % (4 x 3) copy old faces
                edge_faces_def(1,:)
                edge_faces_def(1,:)
                edge_faces_def(2,:)
                edge_faces_def(2,:)
                ];
            new_faces(1,replace_face1(1)) = new_vertex_ID; % (4 x 3) insert new vertex
            new_faces(2,replace_face1(2)) = new_vertex_ID; 
            new_faces(3,replace_face2(1)) = new_vertex_ID; 
            new_faces(4,replace_face2(2)) = new_vertex_ID; 

            faces = [faces;new_faces];           % (#faces x 3)
%    
        
        % Remove 2 old faces
            faces(edge_faces,:) = [];     % (#faces x 3)
%         
        elseif numel(edge_faces)==1
            edge_faces = max(edge_faces);
            new_faces = [            % (4 x 3) copy old faces
                edge_faces_def(1,:)
                edge_faces_def(1,:)
                ];
            new_faces(1,replace_face1(1)) = new_vertex_ID; % (4 x 3) insert new vertex
            new_faces(2,replace_face1(2)) = new_vertex_ID; 
            faces = [faces;new_faces];           % (#faces x 3)

            
        
        % Remove 2 old faces
            faces(edge_faces,:) = [];     % (#faces x 3)

            

        else
            error('More than two faces (%s) on an edge. Check the validity of the triangulation.',mat2str(edge_faces))
        end    

        

end