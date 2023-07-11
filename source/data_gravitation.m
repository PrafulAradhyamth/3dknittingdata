function feature_pts = data_gravitation(V, F, r0, a, b)

% initialize algorithm
%addpath(genpath('..\Extras\geodesic_matlab-master\matlab'))
mesh = geodesic_new_mesh(V, F);     
alg = geodesic_new_algorithm(mesh, 'exact');

% compute gaussian curvature and normalize
k1 = discrete_gaussian_curvature(V, F);
k = (k1 - min(k1))/(max(k1) - min(k1));

% initialize feature points (indices), omit the boundary
bnd = outline(F);
bnd = unique(bnd(:));
interior = setdiff((1:size(V,1))',bnd);
[~, mki] = max(k(interior));
feature_pts = interior(mki);
int_bool = ismember((1:size(V,1))', interior);
% safety break
iters = 0;
% MAIN LOOP
while 1
    iters = iters + 1;
%     fprintf('iteration %d\n', iters)
    % set feature points as sources for the geodesic distance computation
    src = cell(length(feature_pts),1);
    for i=1:length(feature_pts)
        src{i} = geodesic_create_surface_point('vertex',feature_pts(i),V(feature_pts(i),:)); 
    end
    
    % propagate and compute distances
    geodesic_propagate(alg, src,[],5*r0);
    [~, dist] = geodesic_distance_and_source(alg);
    
    % add point with max distance to feature points and check if r < r0
    mask = (dist < 10000) & int_bool;
    mask_idx = find(mask);
    
    [r, mri] = max(dist(mask));
    mr = mask_idx(mri);
    fprintf('iteration %d, max. geodesic distance %f, feat. points len %d\n', iters, r, length(feature_pts))
    feature_pts = cat(1, feature_pts, mr);
%     fprintf('feat. points length %d\n', length(feature_pts))
    if r < r0
        break
    end
    if iters > 150
        warning('Too many iterations! Exiting loop!')
        break
    end
    
    % compute neighbors of every point (within r)
    A = adjacency_list(F);
    
    % compute gravitation for every point
    G = zeros(length(A),1);
    for i = 1:length(A)
        
        if ismember(i, bnd)
            continue;
        end
        
        if k(i) > a
            cx = (k(i)+a)^b;
        else
            cx = k(i);
        end
        G(i) = sum( cx * k(A{i})/r^2 );
    end
    [~, mf] = max(G);
    feature_pts = cat(1, feature_pts, mf);
end