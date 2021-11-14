% create polygons from user-selected tiff
% Note: best to use Mask 1 by default, unless it's a bad mask
% If possible, use same mask for both Cy3 and Cy5
% scripts used: colocalization_v3.m, maskingGUI.m, im2mask.m
% Default masking method: 1 (best for non-noisy images)
% Mask 5 best for early time points/noisy background
% Written by Sophie Skanchy in Version R2019b
% Last updated 9/15/2020

close all

%% >> Input Parameters <<
PixelSize = 0.13; % length of pixel in microns
min_area_um = 0.5; % min droplet area in microns
max_area_um = Inf; % max droplet area in microns

filterCy3_bright = 0; % equalize bright spots from Cy3
filterCy5_bright = 0; % equalize bright spots from Cy5

filterSize = 1;
if filterSize
    min_area = floor(min_area_um / (PixelSize^2)); % convert microns^2 to pixels
    max_area = floor(max_area_um / (PixelSize^2)); % convert microns^2 to pixels
else
    min_area = 0;
    max_area = 1000;
end


%% get tif file
[fname, pathname] = uigetfile('.tif');
cd(pathname);
info = imfinfo(fname);

Cy3_frame = length(info) - 1; % Cy3 frame in TIF
Cy5_frame = length(info); % Cy5 frame in TIF

Cy3 = imread(fname, Cy3_frame, 'Info', info);
Cy5 = imread(fname, Cy5_frame, 'Info', info);
img_size = size(Cy3);

if filterCy3_bright
    Cy3 = medfilt2(adapthisteq(Cy3));
end

if filterCy5_bright
    Cy5 = medfilt2(adapthisteq(Cy5));
end

%% tif to x_ROIs and y_ROIs
% GUI Masking
[Cy3_mask, Cy3_method] = maskingGUI(Cy3);
[Cy5_mask, Cy5_method] = maskingGUI(Cy5);

%% display masks and overlap
disp(fname)
disp(['Used Mask ' num2str(Cy3_method) ' for Cy3, Mask ' num2str(Cy5_method) ' for Cy5'])

% red and green channels
RG = uint8(zeros(img_size(1), img_size(2), 3));
RG(:,:,2) = im2uint8(imadjust(Cy3));
RG(:,:,1) = im2uint8(imadjust(Cy5));

% figure, imshow(Cy3_mask), title('Cy3')
% figure, imshow(Cy5_mask), title('Cy5')
figure, imshow(RG), title('Images overlaid')
saveas(gcf, ['Img_Overlay_' fname(1:4) '.png'])

figure, imshow(imfuse(Cy3_mask, Cy5_mask)), title('Masks overlaid')
saveas(gcf, ['Mask_prefilter_overlay_' fname(1:4) '.png'])

%% label masks
[label_Cy3, numCy3] = bwlabel(Cy3_mask);
[label_Cy5, numCy5] = bwlabel(Cy5_mask);

overlap = Cy3_mask + Cy5_mask > 1;
[~, numOverlap] = bwlabel(overlap);

coloc_fraction = numOverlap / (numCy3 + numCy5 - numOverlap);
disp(['Before filter, coloc fraction: ' num2str(coloc_fraction)])
disp(['Before filter, num Cy3: ' num2str(numCy3)])
disp(['Before filter, num Cy5: ' num2str(numCy5)])
disp(['Before filter, num overlap: ' num2str(numOverlap)])

%% area filter 
    
if filterSize
    % create filter copy of masks
    Cy3_mask_filt = Cy3_mask;
    Cy5_mask_filt = Cy5_mask;
end

Cy3_areas = [];
Cy5_areas = [];

for i = 1:numCy3
    [r, c] = find(label_Cy3 == i);
    Cy3_areas = [Cy3_areas size(r,1)];
    if filterSize && (size(r, 1) < min_area || size(r, 1) > max_area)
        Cy3_mask_filt(r, c) = 0;
    end
end

for i = 1:numCy5
    [r, c] = find(label_Cy5 == i);
    Cy5_areas = [Cy5_areas size(r,1)];
    if filterSize && (size(r, 1) < min_area || size(r, 1) > max_area)
        Cy5_mask_filt(r, c) = 0;
    end
end


if filterSize
    % display filtered masks
    figure, imshow(imfuse(Cy3_mask_filt, Cy5_mask_filt)), title('Masks after area filter overlaid')
    saveas(gcf, ['Mask_filtered_overlay_' fname(1:4) '.png'])
    
    % count number of overlapping droplets
    [~, numCy3] = bwlabel(Cy3_mask_filt);
    [~, numCy5] = bwlabel(Cy5_mask_filt);
    overlap = Cy3_mask_filt + Cy5_mask_filt > 1;
    [~, numOverlap] = bwlabel(overlap);
    coloc_fraction = numOverlap / (numCy3 + numCy5 - numOverlap);
end



%% display results

disp(['Cy3 Droplets: ' num2str(numCy3)])
disp(['Cy5 Droplets: ' num2str(numCy5)])
disp(['Overlapping Droplets: ' num2str(numOverlap)])
disp(['Droplet colocalization fraction: ' num2str(coloc_fraction)])

%% convert areas to microns
Cy3_areas_um = Cy3_areas .* (PixelSize^2);
Cy5_areas_um = Cy5_areas .* (PixelSize^2);

%% histogram of droplet areas before filter
figure, subplot(1,2,1)

% Cy3
Cy3_bins = 0:0.25:ceil(max(Cy3_areas_um));
subplot(1,2,1), histogram(Cy3_areas_um, Cy3_bins)
title('Cy3 Areas')
xlim([0 10])
ylim([0 200])
xlabel('Area (microns^2)')
ylabel('Frequency')
xticks(1:1:10)

% Cy5
Cy5_bins = 0:0.25:ceil(max(Cy5_areas_um));
subplot(1,2,2), histogram(Cy5_areas_um, Cy5_bins)
title('Cy5 Areas')
xlim([0 10])
ylim([0 200])
xlabel('Area (microns^2)')
xticks(1:1:10)

saveas(gcf, ['Area_Distr_' fname(1:4) '.png'])
