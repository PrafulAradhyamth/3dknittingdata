VF = readObj('Zylinder_gebaeugt_45deg_F360.obj');
V = VF.v;
F = VF.f.v;
triangulation(F,V);
trisurf(F,V(:,1), V(:,2), V(:,3));
ptCloud = pointCloud(ans.Points);
pcshow(ptCloud)

% x = ptCloud.Location(:,1);
% y = ptCloud.Location(:,2);
% z = ptCloud.Location(:,3);
% [X, Y] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
% 
% % Interpolate Z values on the grid
% Z = griddata(x, y, z, X, Y, 'linear');
% 
% % Plot the contour
% contour(X, Y, Z);
% 
% % Customize the plot
% title('Contour Plot from Point Cloud');
% xlabel('X');
% ylabel('Y');

