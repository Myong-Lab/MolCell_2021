% Converts .tif of droplets to BW image, counts number of objects
% (droplets), determines area, perimeter, and circularity
% Note: models droplets as ellipses for centroid, major and minor axis, and
% circularity; ellipses plotted with no angle so may not align with droplet
% boundaries, but should have same size
% Counts tend to underestimate number of droplets
% Saves: image-mask comparison (png), image-boundaries comparison (png) 

% parameters: PixelSize = pixel size in um for all images in folder
%             img = image to be analyzed (grayscale)
%             filename = name of .tif, used for saving .png files

% Requires maskings script: im2mask.m, maskingGUI.m and maskingGUI.fig

function [num_droplets, area_stats, circ_stats, img_stats_table] = analyzeDroplets(PixelSize, img, filename)
%% Get parameters for image
% PixelSize = 1.0; %in um

%% convert img to BW mask
% [filename, pathname] = uigetfile('*.tif');
% cd(pathname)
% img = imread(filename);

mask = maskingGUI(img);
img_size = size(img);
saveas(gcf, [filename, '_img2mask', '.png']);

%% regionprops on mask
mask_stats_table = regionprops('table', mask, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Perimeter');
mask_stats_array = table2cell(mask_stats_table);

%% regionprops on mask + grayscale image
%Note: requires BW image of same dimensions as img to use as base
img_stats_table = regionprops('table', mask, img, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'WeightedCentroid', 'MaxIntensity', 'MeanIntensity', 'MinIntensity');
img_stats_array = table2cell(img_stats_table);
img_droplets = length(img_stats_array);
% save img_stats_table

%% label objects in mask
[~, label_droplets] = bwlabel(mask);

%% determine num droplets
if (label_droplets == img_droplets)
    num_droplets = img_droplets;
    disp(['Number of droplets: ' num2str(num_droplets)])
else
    disp('Disagreement between droplet methods. Using count from mask + grayscale image.')
    num_droplets = img_droplets;
    disp(['Number of droplets: ' num2str(num_droplets)])
end

%% arrays for each property
area = zeros(num_droplets, 1);
major_axis_length = zeros(num_droplets, 1);
minor_axis_length = zeros(num_droplets, 1);
perimeter = zeros(num_droplets, 1);
circularity = zeros(num_droplets, 1);
mean_intensity = zeros(num_droplets, 1);
min_intensity = zeros(num_droplets, 1);
max_intensity = zeros(num_droplets, 1);

for i = 1:num_droplets
    area(i) = mask_stats_array{i,1};
    major_axis_length(i) = mask_stats_array{i,3};
    minor_axis_length(i) = mask_stats_array{i,4};
    perimeter(i) = mask_stats_array{i,5};
    circularity(i) = minor_axis_length(i) / major_axis_length(i);
    mean_intensity(i) = img_stats_array{i,6};
    min_intensity(i) = img_stats_array{i,7};
    max_intensity(i) = img_stats_array{i,8};
end

% centroid x and y
centroid_x = zeros(num_droplets, 1);
centroid_y = zeros(num_droplets, 1);
for i = 1:num_droplets
    centroid_x(i) = mask_stats_array{i, 2}(1);
    centroid_y(i) = img_size(1) - mask_stats_array{i, 2}(2);
end

% weighted centroid x and y
w_centroid_x = zeros(num_droplets, 1);
w_centroid_y = zeros(num_droplets, 1);
for i = 1:num_droplets
    w_centroid_x(i) = img_stats_array{i, 5}(1);
    w_centroid_y(i) = img_size(1) - img_stats_array{i, 5}(2);
end

% Convert pixels to um
area = area.*(PixelSize^2);
perimeter = perimeter.*PixelSize;
major_axis_length = major_axis_length.*PixelSize;
minor_axis_length = minor_axis_length.*PixelSize;
centroid_x = centroid_x.*PixelSize;
centroid_y = centroid_y.*PixelSize;
w_centroid_x = w_centroid_x.*PixelSize;
w_centroid_y = w_centroid_y.*PixelSize;

all_props = [area, perimeter, circularity, centroid_x, centroid_y, w_centroid_x, w_centroid_y, major_axis_length, minor_axis_length];

%% plot ellipses representing droplets, compare to plotted ROIs
angle = 0; step = 20;
figure('units','normalized','outerposition',[.1 .1 .9 .7])
title('Compare Input and Output')
subplot(1, 2, 1)
title('ROIs (red) and ellipses (blue)')
for i = 1:num_droplets
    [X,Y] = calculateEllipse(centroid_x(i), centroid_y(i), major_axis_length(i)/2, minor_axis_length(i)/2, angle, step);
    %plot(centroid_x(i), centroid_y(i), 'o')
    %hold on;
    plot(X, Y, '-b')
    hold on;
end
axis([0 img_size(2)*PixelSize 0 img_size(1)*PixelSize])
xlabel('X (um)')
ylabel('Y (um)')
pbaspect([img_size(2)*PixelSize img_size(1)*PixelSize 1])

subplot(1, 2, 2)
imshow(imadjust(img));
title('original tif')
saveas(gcf, [filename, '_ellipses', '.png']);


%% find coords of droplet boundaries
x_ROI = [];
y_ROI = [];
boundaries = bwboundaries(mask);

for i=1:length(boundaries)
    x_ROI=[x_ROI boundaries{i,1}(:,2)'];
    x_ROI(1,end+1)=NaN;
    y_ROI=[y_ROI img_size(1,1)-boundaries{i,1}(:,1)'];
    y_ROI(1,end+1)=NaN;
end

% add ROIs to ellipses plot
subplot(1, 2, 1)
plot(x_ROI*PixelSize, y_ROI*PixelSize, 'r-')


%% calculate stats
% look at: area, count, circularity
% For everything: average, stdev, standard error of mean (sem)
% For count: normalize to count per image
% For area: also include raw area for all droplets

average_area = mean(area);
stdev_area = std(area);
sem_area = stdev_area / sqrt(length(area));
area_stats = [average_area, stdev_area, sem_area];

average_circ = mean(circularity);
stdev_circ = std(circularity);
sem_circ = stdev_circ / sqrt(length(circularity));
circ_stats = [average_circ, stdev_circ, sem_circ];

%% display results
disp(['Average Area: ' num2str(average_area) ' +/- ' num2str(sem_area)])
disp(['Average Circularity: ' num2str(average_circ) ' +/- ' num2str(sem_circ)])
disp(' ')
