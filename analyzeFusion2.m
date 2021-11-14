% Input: .tif video of two droplets fusing
% Output: time constant of relaxation, length scale of fusing droplets
% Parameters: interval, scale, set_ARf, save_mask_vid, gif_speed, remove_zeros

% Uses ellipse model to determine AR
% Uses tau, ARo, ARf as fit parameters
% Required files: im2mask.m, getEllipseAspectRatio.m, masks2gif.m, labelpoints.m
% Note: im2mask.m can be configured with desired masking method

% Written by Sophie Skanchy in version R2019b

close all; clear; clc
%% >> SET PARAMETERS <<
name_abbrev_num = 3; % number of characters from start of filename to include in graph titles
mask_method = 8; % see im2mask.m for method descriptions
invert_img = 1; % set to 1 to invert image (black to white and white to black)
save_plot = 0; % save png of fitted plot
save_mat = 0; % save variables (AR, time, fit) in .mat file

interval = 0.06; % acquisition interval (seconds per frame)
scale = 0.0871; % microns per pixel scale of video


remove_zeros = 0; % removes points where AR = 0 before fitting (default off)

% manual and auto clipping (manual is recommended)
% recommended settings: manual start and end clipping (set autos to zero)
clip_start_auto = 0; % exclude frames until droplets touch; if set to off, fit will not be performed
clip_end_auto = 0; % exclude frames clip_end_length after starting frame
clip_end_auto_length = 6; % number of frames to include (starting at start point)

% mask to animated GIF
save_mask_vid = 1; % save masks as animated gif
gif_speed = 1; % speed multiplier for mask gif

%% get tif file
[fname, pathname] = uigetfile('.tif');
cd(pathname);
info = imfinfo(fname);
vid_length = length(info);
aspect_ratios = zeros(1, vid_length);
ellipse_props = cell(1, vid_length);

name_abbrev = [fname(1:name_abbrev_num) '_v2'];

%% loop through image stack
img_stack = cell(1, vid_length);
masks = cell(1, vid_length);
figure('units','normalized','outerposition',[0 0 1 1])
set(gcf, 'color', 'w')
for i = 1:vid_length
    img = imread(fname, i, 'Info', info);
    if invert_img
        img = imcomplement(img);
    end
    img_stack{i} = img;
    img_size = size(img);
    mask = im2mask(img, mask_method);
    masks{i} = mask;
    %figure
    imshow(imfuse(imadjust(img), mask))
    [aspect_ratios(i), ellipse_props{i}] = getEllipseAspectRatio(mask);
    clc; disp(['Analyzing Frame: ' num2str(i)])
end

%% save GIF of masks
if save_mask_vid == 1
    mask_gif_name = [fname(1:end-4) '_masks.gif'];
    masks2gif(masks, img_stack, interval, mask_gif_name, gif_speed);
end

%% manually select start and end points
if clip_start_auto == 0 || clip_end_auto == 0
    disp('Plotting aspect ratios...')
    frame = 1:length(aspect_ratios);
    figure, scatter(frame, aspect_ratios), title('Aspect Ratio')
    xlabel('frame'); ylabel('Aspect Ratio');
    disp('Labeling data points...')
    labels = cell(1, length(frame));
    for i = 1:length(frame)
        labels{i} = num2str(frame(i));
    end
    labelpoints(frame, aspect_ratios, labels); % label data points
    if clip_start_auto == 0
        start_frame = input('Choose start frame (#): ');
        while start_frame < 1 || start_frame > length(aspect_ratios)
            start_frame = input('Invalid start frame. Reenter start frame: ');
        end
    end
    if clip_end_auto == 0
        end_frame = input('Choose end frame (#): ');
        while end_frame < 2 || end_frame > length(aspect_ratios)
            end_frame = input('Invalid start frame. Reenter end frame: ');
        end
    end
end

%% automatically clip start and determine start point
if clip_start_auto
    for i = 1:vid_length
        if aspect_ratios(i) ~= 0 && aspect_ratios(i+1) ~= 0
            ARo = aspect_ratios(i);
            ind_ARo = i;
            break
        end
    end
    start_frame = ind_ARo;
    aspect_ratios = aspect_ratios(ind_ARo:vid_length);
elseif ~clip_start_auto && start_frame > 1
    ARo = aspect_ratios(start_frame);
    ind_ARo = start_frame;
    aspect_ratios = aspect_ratios(start_frame:vid_length); % clip aspect_ratios
else
    start_frame = 1;
    ARo = aspect_ratios(1);
    ind_ARo = 1;
end
time = 0:length(aspect_ratios) - 1;
time = time * interval;


%% clip end
if clip_end_auto
    end_frame = start_frame + clip_end_auto_length;
    aspect_ratios = aspect_ratios(1:clip_end_auto_length);
    time = time(1:clip_end_auto_length);
elseif ~clip_end_auto && end_frame < vid_length
    aspect_ratios = aspect_ratios(1:end_frame - start_frame + 1);
    time = time(1:end_frame - start_frame + 1);
end


%% plot and fit aspect ratios over time
% remove zeros
if remove_zeros == 1
    nonzero_inds = aspect_ratios ~= 0;
    aspect_ratios = aspect_ratios(nonzero_inds);
    time = time(nonzero_inds);
end

% plot and fit
figure, scatter(time, aspect_ratios, 'k.')
fit_func = 'ARf+(ARo-ARf)*exp(-time/tau)';
fit_type = fittype(fit_func, 'independent',{'time'}, 'coefficients',{'tau', 'ARf', 'ARo'});
[fitted_curve, gof] = fit(time', aspect_ratios', fit_type, 'StartPoint', [0.1, 2, 1]);
hold on, plot(fitted_curve)
xlabel('Time (s)'), ylabel('Aspect Ratios')
legend('data', 'fitted curve')
title([name_abbrev '.tif Fusion: Aspect Ratio vs. Time'])
%save figure, include tif # in filename
ARf = fitted_curve.ARf;
conf_ints = confint(fitted_curve);
ARf_range = conf_ints(:,2);
delta_ARf = abs(ARf - ARf_range(1));

ARo = fitted_curve.ARo;
conf_ints = confint(fitted_curve);
ARo_range = conf_ints(:,3);
delta_ARo = abs(ARo - ARo_range(1));

% display tau with uncertainty and r^2
tau = fitted_curve.tau;
intervals = confint(fitted_curve);
tau_range = intervals(:,1);
delta_tau = abs(tau - tau_range(1));
annotation('textbox',[.59 .5 .31 .28],'String',sprintf('tau = %.4f +/- %.3f s', tau, delta_tau),'EdgeColor','k')
annotation('textbox',[.6 .65 .31 .1],'String',sprintf('r^2 = %.5f', gof.rsquare),'EdgeColor','none')
annotation('textbox',[.6 .55 .31 .1],'String',sprintf('ARo = %.4f +/- %.4f', ARo, delta_ARo),'EdgeColor','none')
annotation('textbox',[.6 .52 .31 .1],'String',sprintf('ARf = %.4f +/- %.4f', ARf, delta_ARf),'EdgeColor','none')
annotation('textbox',[.6 .48 .31 .1],'String',sprintf('start frame: %d', start_frame),'EdgeColor','none')
annotation('textbox',[.6 .45 .31 .1],'String',sprintf('end frame: %d', end_frame),'EdgeColor','none')

%% calculate length scale
long = ellipse_props{start_frame}(1).MajorAxisLength*scale;
short = ellipse_props{start_frame}(1).MinorAxisLength*scale;
length_scale = sqrt((long - short)*short);
annotation('textbox',[.6 .63 .31 .07],'String',sprintf('Length Scale: %.2f um', length_scale),'EdgeColor','none')

%% display results
disp(' ')
disp(['Fusion analysis for ' fname])
disp(['Start frame: ' num2str(start_frame)])
disp(['End frame: ' num2str(end_frame)])
disp(['ARo: ' num2str(ARo)])
disp(['ARf: ' num2str(ARf)])
disp(['Time constant: ' num2str(tau) ' +/- ' num2str(delta_tau)])
disp(['r square: ' num2str(gof.rsquare)])
disp(['Length scale: ' num2str(length_scale)])

%% save figure
if save_plot
    saveas(gcf, [name_abbrev '_plot.png'])
end

if save_mat
    save([name_abbrev '.mat'], 'mask_method', 'invert_img', 'start_frame', ...
        'end_frame', 'time', 'aspect_ratios', 'fitted_curve', 'gof', 'tau', ...
        'delta_tau', 'ARo', 'delta_ARo', 'ARf', 'delta_ARf', 'length_scale')
end

