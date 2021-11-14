% Determine aspect ratio given BW mask of fused droplets and returns ellipse properties.
% If more than one droplet detected, selects object with index 1.
% Method: model fused droplets as an ellipse, AR = major axis / minor axis
% Written by Sophie Skanchy in R2019b

% @param BW binary image of droplets fusing

function [aspect_ratio, ellipse_props] = getEllipseAspectRatio(BW)
%% body
% regionprops to model as ellipse
multiple_droplets = 0;
no_droplets = 0;
ellipse_props = regionprops('struct', BW, 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
droplet_choice = 1; % if multiple droplets, uses this default index to access droplet

% calculate aspect ratio as major axis/minor axis
if length(ellipse_props) > 1
     disp(['Multiple droplets detected, AR based on droplet #' num2str(droplet_choice)])
     multiple_droplets = 1;
     
     % find largest blob instead of using default index
     maj = zeros(1,length(ellipse_props));
     for i = 1:length(ellipse_props)
        maj(i) = ellipse_props(i).MajorAxisLength;
     end
     droplet_choice = find(maj == max(maj));
     
     aspect_ratio = ellipse_props(droplet_choice).MajorAxisLength / ellipse_props(droplet_choice).MinorAxisLength;;
elseif length(ellipse_props) < 1
    disp('No droplets detected')
    no_droplets = 1;
    aspect_ratio = 0;
else
    aspect_ratio = ellipse_props.MajorAxisLength / ellipse_props.MinorAxisLength;
end

%% plot droplet boundary and ellipse
% BW to boundary coords
% if multiple_droplets == 0 && no_droplets == 0
%     BW = flip(BW);
%     x_ROI = [];
%     y_ROI = [];
%     boundaries = bwboundaries(BW);
%     mask_size = size(BW);
%     for j=1:length(boundaries)
%         x_ROI=[x_ROI boundaries{j,1}(:,2)'];
%         x_ROI(1,end+1)=NaN;
%         y_ROI=[y_ROI mask_size(1,1)-boundaries{j,1}(:,1)'];
%         y_ROI(1,end+1)=NaN;
%     end
    
    % plot
%     angle = ellipse_props.Orientation; step = 20;
%     [eX,eY] = calculateEllipse(ellipse_props.Centroid(1), ellipse_props.Centroid(2), ellipse_props.MajorAxisLength/2, ellipse_props.MinorAxisLength/2, angle, step);
%     scatter(ellipse_props.Centroid(1), ellipse_props.Centroid(2), 'bo')
%     hold on; plot(eX, eY, '-b')
%     hold on; plot(x_ROI, y_ROI, 'r')
%     title('Droplet Boundary and Ellipse Model')
end

