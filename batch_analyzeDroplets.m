% Batches droplet analysis for all .tif files in user-selected directory.
% For each .tif, saves table of values from regionprops analysis, number of
% droplets, and area and circularity statistics. File names of saved data
% based on .tif filename.

close all; clc

%% Get folder with .tif files
Directory = uigetdir();
cd(Directory);
filenames = dir([Directory, '\*.tif']);
number_of_files = length(filenames);
PixelSize = input('Size of pixel in microns: '); % for all images in folder
out_name = input('Output .txt filename: ', 's');
while ~strcmp(out_name(end-3:end), '.txt') % check for valid out_name
    out_name = input('Please input valid .txt filename: ', 's');
end
disp(' ')

%% Create .txt to save data for each image
% fields are comma-separated
%out_name = 'droplet_data.txt';
out_filedir = Directory;

try %try-catch block to close file in case of error
    fileID = fopen(out_name,'w');
    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s\n', 'Filename', 'Number of Droplets', 'Mean Area', 'Area STDEV', 'Area SEM', 'Mean Circularity', 'Circularity STDEV', 'Circularity SEM');
    
    %% loop through .tif files in folder
    droplet_counts = zeros(1, number_of_files);
    all_area = [];
    all_circ = [];
    
    for fileNum = 1:number_of_files
        filename = filenames(fileNum).name;
        img = imread(filename);
        disp(['Analyzing ', filename]);
        [num_droplets, area_stats, circ_stats, img_stats_table] = analyzeDroplets(PixelSize, img, filename(1:length(filename)-4));
        droplet_counts(fileNum) = num_droplets;
        all_area = horzcat(all_area, area_stats(1));
        all_circ = horzcat(all_circ, circ_stats(1));
        
        %save desired info (tables as separate files, .txt with appended num_droplets, area_stats, circ_stats)
        writetable(img_stats_table, [filename(1:length(filename)-4), '_stats_table'])
        fprintf(fileID,'%s,%d,%f,%f,%f,%f,%f,%f\n', filename, num_droplets, area_stats(1), area_stats(2), area_stats(3), circ_stats(1), circ_stats(2), circ_stats(3));
    end
    
    fileID = fclose(fileID);
    
    %% Calculate stats for droplet_counts, all_area, all_circ
    % droplet counts
    average_droplets = mean(droplet_counts);
    stdev_droplets = std(droplet_counts);
    sem_droplets = stdev_droplets / sqrt(length(droplet_counts));
    disp('Droplet Count stats ')
    disp(['Average: ' num2str(average_droplets) ', STDEV: ' num2str(stdev_droplets) ', SEM: ' num2str(sem_droplets)])
    
    % area (note: all_area is list of mean areas for each image)
    average_area = mean(all_area);
    stdev_area = std(all_area);
    sem_area = stdev_area / sqrt(length(all_area));
    disp(' ')
    disp('Area stats ')
    disp(['Average: ' num2str(average_area) ', STDEV: ' num2str(stdev_area) ', SEM: ' num2str(sem_area)])
    
    % circularity (note: all_circ is list of mean circularity for each image)
    average_circ = mean(all_circ);
    stdev_circ = std(all_circ);
    sem_circ = stdev_circ / sqrt(length(all_circ));
    disp(' ')
    disp('Circularity stats ')
    disp(['Average: ' num2str(average_circ) ', STDEV: ' num2str(stdev_circ) ', SEM: ' num2str(sem_circ)])
    
catch ME % catch errors so .txt file is always closed
    if fileID ~= 0
        fclose(fileID);
    end
    disp(['Error encountered at ' num2str(ME.stack.line) ': ' ME.identifier])
    disp(ME.message)
    disp('Closed output .txt file')
end
