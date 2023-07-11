function [level_time] = pc_rotate(points)
%points are the level time 
% load('matlab.mat')
% points = level_time{25,1};
% plot3(points(1,:),points(2,:),points(3,:))
% figure;
%taking any three points to find the angle of the plane that these points sit
% with respect to XY plane, (that can be modified)
P1 = points(:,1);
P2 = points(:,100);
P3 = points(:,55);
% Calculate the normal vector of the circle's plane
v1 = P2 - P1;
v2 = P3 - P1;
normal = cross(v1, v2);
normal = normal / norm(normal);
% Define the normal vector of the XY plane
normalXY = [0, 0, 1]; % Assuming the XY plane is the standard plane
% Calculate the dot product between the two normal vectors
dotProduct = dot(normal, normalXY);
% Calculate the angle between the two planes in radians
angleRadXY = acos(dotProduct);
% Convert the angle to degrees
angleDegXY = rad2deg(angleRadXY);

%
normalXZ = [0, 1, 0]; % Assuming the XY plane is the standard plane
% Calculate the dot product between the two normal vectors
dotProduct = dot(normal, normalXZ);
% Calculate the angle between the two planes in radians
angleRadXZ = acos(dotProduct);
% Convert the angle to degrees
angleDegXZ = rad2deg(angleRadXZ);

%
normalYZ = [1, 0, 0]; % Assuming the XY plane is the standard plane
% Calculate the dot product between the two normal vectors
dotProduct = dot(normal, normalYZ);
% Calculate the angle between the two planes in radians
angleRadYZ = acos(dotProduct);
% Convert the angle to degrees
angleDegYZ = rad2deg(angleRadYZ);

%if condition based on the iso lines and the straightness of the 

%rotation of point cloud
center1 = mean(points');
rotationAngles = [0,-angleDegXZ,0];
translation = [0 0 0];
tform = rigidtform3d(rotationAngles,translation);
ptCloudOut = pctransform(pointCloud(points'),tform);
center2 = mean(pointCloud(points').Location);

%rotation of point cloud
rotationAngles = [0,0, 0];
translation = [center2(1)-center1(1), center2(2)-center1(2), center2(3)-center1(3)];
tform = rigidtform3d(rotationAngles,translation);
ptCloudOut = pctransform(pointCloud(points'),tform);

% pcshow(pointCloud(points'))
% figure;
% pcshow(ptCloudOut)
% figure;
%ptCloudOut.Location  
level_time = ptCloudOut.Location ;
%plot3(level_time(:,1),level_time(:,2),level_time(:,3))

end