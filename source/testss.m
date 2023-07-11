% Example point cloud data
x = rand(100, 1) * 10;
y = rand(100, 1) * 10;
z = sin(x) + cos(y);

% Construct the point cloud matrix
point_cloud = [x, y, z];
pcshow(point_cloud)
% Compute the Reeb graph
r = rips(point_cloud, 2, 2); % 2 is the max dimension and 2 is the max distance

% Plot the Reeb graph
plot(r);

% Customize the plot
title('Reeb Graph');
xlabel('Node Index');
ylabel('Persistence');

% Add labels to the nodes
for i = 1:numel(r.nodes)
    text(r.nodes(i).index, r.nodes(i).persistence, num2str(i));
end
