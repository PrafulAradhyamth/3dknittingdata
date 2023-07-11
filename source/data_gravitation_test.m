% test data gravitation
clear; clc; close all
fclose('all');
set_interpreter_paths
addpath(genpath('..\Extras\geodesic_matlab-master\matlab'))

global geodesic_library;                
geodesic_library = 'geodesic_release';

pathstr = '..\Data\';
filenames = {'Last_High_Shaft.obj', 'Last_Low_Shaft.obj',...
        'TEST2.obj','Sitzkissen_surface_remesh.obj',...
        'Sitzschale 2teilig_1_Lehne_surface_remesh.obj',...
        'Sitzschale 2teilig_2_Sitz_surface_remesh.obj'};

for i=1:length(filenames)
    
    [V,F] = readOBJ(strcat([pathstr,filenames{i}]));
    dx = max(V(:,1))-min(V(:,1));
    dy = max(V(:,2))-min(V(:,2));
    dz = max(V(:,3))-min(V(:,3));
    diam = sqrt(dx^2+dy^2+dz^2);
    fprintf('Examining file %s, mesh diameter %f\n',filenames{i}, diam)
    params = [diam/2, 2,4; diam/3, 2,4; diam/4, 2,4; diam/5, 2,4];
    figure;
    for j=1:4
        feature_pts = data_gravitation(V, F, params(j,1), params(j,2), params(j,3));        
        subplot(2,2,j)
        trisurf(triangulation(F,V))%,'facecolor','none','edgecolor',[0.7,0.7,0.7]); 
        axis equal
        hold on
        scatter3(V(feature_pts,1),V(feature_pts,2),V(feature_pts,3),'ro','filled')
        hold off
    end
end

% a,b, have hardly any influence 