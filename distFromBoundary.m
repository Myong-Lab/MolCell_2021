function [minDist, indMinDist] = distFromBoundary(xCoord, yCoord, x_ROI, y_ROI)
%% Description
% Given x_ROI (x coords of droplet boundaries) and y_ROI (y coords of
% droplet boundaries), calculates the distance between a point given by
% xCoord and yCoord and the closest point on the boundaries. Returns this
% distance, and the index in x_ROI and y_ROI of the closest boundary point.

%% Body
% calculate x dist and y dist of point from each point in boundary
dx = x_ROI - xCoord;
dy = y_ROI - yCoord;

% dist = sqrt(x^2 + y^2), from each point in x_ROI, y_ROI
dist = sqrt(dx.^2 + dy.^2);
[minDist, indMinDist] = min(dist);

% % Test:
% figure; plot(x_ROI, y_ROI);
% hold on; scatter(xCoord, yCoord, 'o');
% hold on; scatter(x_ROI(indMinDist), y_ROI(indMinDist), 'o');

%Note: x coord of the closest point on boundary is x_ROI(indMinDist),
%      y coord is y_ROI(indMinDist)