function [] = visODF(datapoints, directions)
%%3D visualization code for points on specified by directions, and 
%%radius weighed by the values in datapoints.
%%Size of datapoints should be same as # of rows in directions.
%%# of columns in directions should strictly be 3 (for 3D)

fcs = convhulln(directions);

points(:,1) = datapoints.*directions(:,1);
points(:,2) = datapoints.*directions(:,2);
points(:,3) = datapoints.*directions(:,3);

patch('Vertices',[points(:,1) points(:,2) points(:,3)],'Faces',fcs,'FaceVertexCData',datapoints,...
    'FaceColor','interp','EdgeColor','interp','BackFaceLighting','lit',...
    'FaceLighting','phong');
grid on;
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');
view([1 1 1]);
colorbar;
