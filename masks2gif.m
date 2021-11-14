% Given cell array of BW masks and time interval, creates and saves animated GIF
% Parameters: BWarray is cell array of BW images
%             interval is acquisition time (time per frame)
%             filename is string, name to save gif with, should end with .gif!
%             speed is multiplier for video speed (e.g. speed = 2 is 2x speed)
% returns 0 if unsuccessful, 1 if successful
% Written by Sophie Skanchy in R2019b

function success = masks2gif(BWcells, img_stack, interval, filename, speed)
vid_length = length(BWcells);
fig = figure;
imcells = cell(1, vid_length);

for i = 1:vid_length
    clc; disp(['Saving to GIF: Frame ' num2str(i)])
    imshow(imfuse(imadjust(img_stack{i}),BWcells{i}))
    frame = getframe(fig);
    imcells{i} = frame2im(frame);
end

for idx = 1:vid_length
    [A,map] = rgb2ind(imcells{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',interval/speed);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',interval/speed);
    end
end
success = 1;