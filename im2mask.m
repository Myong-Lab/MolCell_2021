% converts grayscale image to BW mask
% for use with maskingGUI.m, analyzeDroplets.m 
% FINAL VERSION for DROPLET ANALYSIS as of 3/10/2020

% note: if masking methods are added, must update maskingGUI.fig and maskVideoGUI.fig dropdowns
% Mask methods are optimized for fluorescence images unless otherwise stated

% Masks list:
% Mask 1: Amir's Method
% Mask 2: Sobel Edge
% Mask 3: Canny + log Edge
% Mask 4: log Edge
% Mask 5: Median Filtering and Graythresh
% Mask 6: Majority of masks 1-5, DEFAULT MASK



function mask = im2mask(I, method)
% img is original grayscale image
% method is chosen masking method (1-4)

num_masks = 6;
mask_size = size(I);
mask = cell(1, num_masks);

%% Mask 1: Amir's Method
level=ceil(graythresh(I)*100)/100;
BW=im2bw(I,level);

% cleaning the BW image
%BW = bwmorph(BW, 'clean');
BW2 = imfill(BW, 'holes');
BW3 = bwmorph(BW2,'majority', 2);
BW4 = bwareaopen(BW3, 5);
BW5 = imclearborder(BW4);
BW6 = bwmorph(BW5,'thicken',1);

mask{1} = BW6;


%% Mask 2: Edge Detection - Sobel
edgeBW = edge(I, 'Sobel');
%close gaps
radius = 2; num = 4;
se = strel('disk', radius, num);
edgeBWc = imclose(edgeBW,se);

%fill interior
edgeBWfill = imfill(edgeBWc,'holes');

BW2 = imfill(edgeBWfill, 'holes');
BW3 = bwmorph(BW2,'majority', 2);
BW4 = bwareaopen(BW3, 5);
BW5 = imclearborder(BW4);
BW6 = bwmorph(BW5,'thicken',1);

mask{2} = BW6;


%% Mask 3: Edge Detection - Canny + log, no strel
% edge detection
edgeCanny = edge(I, 'Canny');
edgeLog = edge(I, 'log');

% fill and remove noise
fillCanny = imfill(edgeCanny, 'holes');
majCanny = bwmorph(fillCanny, 'majority');

fillLog = imfill(edgeLog, 'holes');
majLog = bwmorph(fillLog, 'majority');

% sum majCanny and majLog, remove noise
majCannyLog = (majCanny + majLog) > 0;
majCannyLog = bwmorph(majCannyLog,'majority', 3);

mask{3} = majCannyLog;


%% Mask 4: Edge Detection - log, no strel, minimal cleaning
edgeLog = edge(I, 'log');
fillLog = imfill(edgeLog, 'holes');
majLog = bwmorph(fillLog, 'majority', 2);

mask{4} = majLog;


%% Mask 5: Median Filtering and Graythresh
clear = imclearborder(I);
equalize = adapthisteq(clear);
filter = medfilt2(equalize);
level=ceil(graythresh(filter)*110)/100;
BWfilter = im2bw(filter, level);
BWclean = bwmorph(BWfilter,'majority',2);

mask{5} = BWclean;


%% Mask 6: Majority of masks 1-5
sum_masks = mask{1} + mask{2} + mask{3} + mask{4} + mask{5};
average = sum_masks > 2; % at least 3 of 5 masks agree
clean = bwmorph(average,'majority',2);

mask{6} = clean;


%% select mask
mask = mask{method};