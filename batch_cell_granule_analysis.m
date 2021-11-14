% cell image granule analysis
% Written by: Sophie Skanchy
% Last updated: 1/4/2021
% Written in Matlab R2019b
% scripts needed: distFromBoundary.m, image processing toolbox
% Outputs: txt file for batch, txt files for each image with granules list,
% .mat with data and masks, .png with masks and nearest nucs
% all distance outputs are in microns (or um^2 for areas)

close all;
clc;
clear;

%% >> Input Parameters <<
PixelSize = 0.1984666; % microns per pixel

% filters
min_area_um = 200*PixelSize^2; % min nucleus area in microns^2
min_overlap_area_um = 1000*PixelSize^2; % min overlap between nuclei and FUS in microns^2

max_gran_dist_um = 20; % max distance of granule from nucleus in microns
max_gran_area_um = 10; % max area of granules in microns^2

% masking
FUS_mask_method = 3; % options: 1 thru 5, 4 and 5 old/not optimal
gran_mask_method = 3; % options: 1 thru 4, 4 not optimal

% misc
watershed_slicing = 0; % 0 for no, 1 for yes; for nuclei mask
remove_nuc_gran = 1; % 0 for no, 1 for yes, only removes if entire granule overlaps nucleus

out_name = 'batch_stats.txt';
while ~strcmp(out_name(end-3:end), '.txt') % check for valid out_name
    out_name = input('Please input valid .txt filename: ', 's');
end
disp(' ')

% do not change
min_area = min_area_um/PixelSize^2; % minimum nucleus area in pixels
min_overlap_area = min_overlap_area_um/PixelSize^2; % minimum overlap area in pixels
max_gran_dist = max_gran_dist_um/PixelSize; % max distance of granule from nucleus in pixels
max_gran_area = max_gran_area_um/PixelSize^2; % max area of granules in pixels

%% get directory
PathName = uigetdir;
cd(PathName);

filename=dir([PathName, '\*.tif']);
number_of_files=length(filename);

%% Create batch .txt
try %try-catch block to close file in case of error
    fileID = fopen(out_name,'w');
    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'Filename', 'Num_Granules', ...
        'Total_Area', 'Mean_Area', 'Area_STDEV', 'Area_SEM', 'Num_Nuclei', ...
        'Num_GranulesNF', 'TotalAreaNF', 'Mean_AreaNF', 'Area_StDevNF', 'Area_SEMNF');
    % NF means no granule filters
    
    %% loop through files
    for f = 1:number_of_files
        %% import image frames
        file = filename(f).name;
        
        img_info = imfinfo(file);
        DAPI_o = imread(file, 2, 'Info', img_info); % nuclei
        FUS_o = imread(file, 1, 'Info', img_info); % FUS
        
        FUS = imadjust(FUS_o); % Green
        DAPI = imadjust(DAPI_o); % Blue
        
        origRGB = cat(3, DAPI_o, zeros(size(FUS_o)), FUS_o);
        figure, imshow(origRGB), title('orig rgb')
        %saveas(gcf, [file(1:3), '_orig', '.png']);
        
        disp(['Analyzing ' file])
        
        %% mask DAPI
        % imbinarize
        sens = 0.6;
        BWa = imbinarize(imgaussfilt(DAPI_o), 'adaptive', 'Sensitivity', sens); % adaptive
        %figure, imshow(imfuse(orig, BWa)), title('adaptive imbinarize')
        
        % clean binarized image
        BW2 = imfill(BWa, 'holes');
        BW3 = bwmorph(BW2,'majority', 5);
        BW4 = bwareaopen(BW3, 5);
        BW5 = imclearborder(BW4);
        clean_binarize = bwmorph(BW5,'thicken',1);
        
        %figure, imshow(clean_binarize), title('cleaned adaptive imbinarize')
        
        % graythresh
        I = DAPI;
        level=ceil(graythresh(I)*100)/100;
        orig_thresh = im2bw(I,level);
        
        %figure, imshow(orig_thresh), title('orig > graythresh')
        
        % combine graythresh and imbinarize
        combineBW = orig_thresh + clean_binarize > 1;
        combineBW = imfill(combineBW, 'holes');
        combineBW = bwmorph(combineBW,'majority', 5);
        DAPI_mask = combineBW;
        %figure, imshow(combineBW), title('graythresh + adapt imbinarize')
        
        % watershed segmentation
        if watershed_slicing
            bw = combineBW;
            D = -bwdist(~bw);
            watershed_mask = imextendedmin(D,2);
            D2 = imimposemin(D,watershed_mask);
            Ld2 = watershed(D2);
            bw2 = bw;
            bw2(Ld2 == 0) = 0;
            DAPI_mask = bw2;
        end
        
        %figure, imshow(DAPI_mask), title('imbinarize+graythresh > watershed')
        figure, imshow(imfuse(DAPI, DAPI_mask)), title('DAPI mask')
        
        %% mask FUS
        if FUS_mask_method == 1 % Method 1: low pass and equalize > graythresh
            I_eq = adapthisteq(imgaussfilt(FUS_o));
            bw = im2bw(I_eq, graythresh(I_eq));
            bw2 = imfill(bw,'holes');
            bw3 = imopen(bw2, ones(5,5));
            FUS_mask = bwareaopen(bw3, 40);
            %figure, imshow(imfuse(FUS, FUS_mask)), title('FUS mask')
            
        elseif FUS_mask_method == 2 % Method 2: equalize and high pass > graythresh
            K = adapthisteq(FUS_o) - imgaussfilt(FUS_o, 2.5); % high pass
            bw = im2bw(K, graythresh(K));
            bw2 = imfill(bw,'holes');
            bw3 = imopen(bw2, ones(5,5));
            FUS_mask = bwareaopen(bw3, 40);
            %figure, imshow(FUS_mask)
            
        elseif FUS_mask_method == 3 % Method 3: union of methods 1 and 2
            % Method 1
            I_eq = adapthisteq(imgaussfilt(FUS_o));
            bw = im2bw(I_eq, graythresh(I_eq));
            bw2 = imfill(bw,'holes');
            bw3_1 = imopen(bw2, ones(5,5));
            
            % Method 2
            K = adapthisteq(FUS_o) - imgaussfilt(FUS_o, 2.5); % high pass
            bw = im2bw(K, graythresh(K));
            bw2 = imfill(bw,'holes');
            bw3_2 = imopen(bw2, ones(5,5));
            
            FUS_mask = bw3_1 + bw3_2 > 0;
            
            % METHODS BELOW ARE OLDER/LESS OPTIMAL
        elseif FUS_mask_method == 4 % Method 4: graythresh
            % graythresh
            I = FUS;
            threshold_adjust = 100; % percent multiplier for intensity threshold
            level=ceil(graythresh(I)*threshold_adjust)/100;
            orig_thresh = im2bw(I,level);
            
            % clean
            FUS_mask = bwmorph(orig_thresh,'majority', 2);
            
        elseif FUS_mask_method == 5 % Method 5: imbinarize + graythresh
            % imbinarize
            sens = 0.6;
            BWa = imbinarize(imgaussfilt(FUS_o), 'adaptive', 'Sensitivity', sens); % adaptive
            
            % clean binarized image
            BW2 = imfill(BWa, 'holes');
            BW3 = bwmorph(BW2,'majority', 2);
            BW4 = bwareaopen(BW3, 5);
            BW5 = imclearborder(BW4);
            clean_binarize = bwmorph(BW5,'thicken',1);
            
            % graythresh
            I = FUS;
            threshold_adjust = 100; % percent multiplier for intensity threshold
            level=ceil(graythresh(I)*threshold_adjust)/100;
            orig_thresh = im2bw(I,level);
            
            % clean
            clean_thresh = bwmorph(orig_thresh,'majority', 2);
            
            % combine graythresh and imbinarize
            FUS_mask = clean_thresh + clean_binarize > 1;
        end
        
        
        % Ignore:
        % cell_outline = bwperim(FUS_mask);
        % Segout = I;
        % Segout(cell_outline) = 255;
        % figure, imshow(Segout)
        % title('Outlined Cell')
        
        %% mask FUS granules
        
        if gran_mask_method == 1 % Method 1: imadjust then high pass
            % adjust contrast
            low_in = 0.15; high_in = 1;
            low_out = 0; high_out = 1;
            gamma = 1.3;
            
            J = imadjust(FUS_o,[low_in high_in],[low_out high_out], gamma);
            
            % subtract background (high pass filter)
            sigma = 2.5; % for low pass/Gaussian filter
            BG = imadjust(imgaussfilt(FUS_o, sigma), [0 0.8], [0 0.2]); % adjusted low pass
            K = J - BG; % high pass = image - low pass
            
            % graythresh
            level = ceil(graythresh(K)*100)/100;
            BW = im2bw(K, level);
            clean_granules = bwmorph(BW,'majority', 2);
            
        elseif gran_mask_method == 2 % Method 2: imadjust after high pass (less sensitive)
            low_pass = imgaussfilt(FUS_o, sigma);
            high_pass = FUS_o - low_pass;
            K = imadjust(high_pass,[.001 .999], [0 1], gamma);
            
            % graythresh
            level = ceil(graythresh(K)*100)/100;
            BW = im2bw(K, level);
            clean_granules = bwmorph(BW,'majority', 3);
        elseif gran_mask_method == 3 % Method 3: normalize high pass
            norm = uint8(255*mat2gray(FUS_o));
            adjust = imadjust(norm, [0.2 1], [0 1], 1.3);
            blur = imgaussfilt(norm, 20);
            BW = imbinarize(adjust - blur);
            clean_granules = bwmorph(BW, 'majority', 1);
            
        elseif gran_mask_method == 4 % Method 4: subtract mean high pass, backup for unstressed
            norm = uint8(255*mat2gray(FUS_o));
            adjust = imadjust(norm, [0.2 1], [0 1], 1.3);
            thresh = round(mean(adjust,'all'));
            BG = (adjust > thresh).*thresh;
            blur = imgaussfilt(norm, 5);
            
            filter = adjust - blur - uint8(BG);
            filter(filter < 0) = 0;
            BW = imbinarize(filter);
            clean_granules = bwmorph(BW, 'majority', 1);
            
        elseif gran_mask_method == 5 % Method 5: extended maxima (old, less optimal)
            granules = imextendedmax(imgaussfilt(FUS), 1);
            BW = bwmorph(granules, 'majority', 5);
            BW2 = imfill(BW, 'holes');
            BW3 = bwareaopen(BW2, 5);
            clean_granules = imclearborder(BW3);
        end
        
        %figure, imshow(clean_granules)
        figure, imshow(imfuse(FUS, clean_granules)), title('granules mask')
        
        
        %% remove nuclei not overlapping FUS mask
        % count and label nuclei
        [DAPI_label, num_nuclei] = bwlabel(DAPI_mask);
        
        % determine nuclei area and filter by min_area
        nuclei_areas = zeros(1, num_nuclei);
        DAPI_filter = DAPI_label;
        
        for i = 1:num_nuclei
            [r, c] = find(DAPI_label == i);
            nuclei_areas(i) = size(r, 1);
            if size(r, 1) < min_area
                DAPI_filter(r, c) = 0; % filter nuclei by area
            end
        end
        
        % find nuclei numbers overlapping FUS, filter overlap area
        overlap = DAPI_mask + FUS_mask > 1;
        [overlap_label, num_overlap] = bwlabel(overlap);
        overlap_areas = zeros(1, num_overlap);
        
        FUS_filter = FUS_mask;
        for i = 1:num_overlap
            [r, c] = find(overlap_label == i);
            overlap_areas(i) = size(r, 1);
            if size(r, 1) < min_overlap_area
                FUS_filter(r, c) = 0;
            end
        end
        
        
        % get indices of nuclei in FUS mask
        ind_nuclei = unique(DAPI_filter(FUS_filter == 1));
        if ind_nuclei(1) == 0
            ind_nuclei = ind_nuclei(2:end); % remove zero index
        end
        
        % apply overlap filter to DAPI mask
        for r = 1:size(DAPI_label, 1)
            for c = 1:size(DAPI_label, 2)
                if ~ismember(DAPI_label(r, c), ind_nuclei)
                    DAPI_filter(r, c) = 0;
                end
            end
        end
        
        % display number of nuclei after filter
        [~, num_overlap] = bwlabel(DAPI_filter > 0);
        disp(['Nuclei in FUS mask: ' num2str(num_overlap)])
        %         figure, imshow(imfuse(FUS_mask, DAPI_filter > 0))
        %         title('Nuclei in FUS only')
        
        % regionprops area filtering of nuclei
        % DAPI_overlap_props = regionprops('table', overlap, 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Centroid', 'Orientation');
        % ind_keep_overlap = DAPI_overlap_props.Area >= min_area;
        
        %% remove cell bodies in FUS mask not overlapping nuclei
        % label FUS mask
        FUS_label = bwlabel(FUS_mask);
        
        % get indices of FUS in DAPI mask
        ind_FUS = unique(FUS_label(DAPI_filter >= 1));
        if ind_FUS(1) == 0
            ind_FUS = ind_FUS(2:end); % remove zero index
        end
        
        % apply overlap filter to FUS mask
        FUS_filter = FUS_mask;
        for r = 1:size(FUS_label, 1)
            for c = 1:size(FUS_label, 2)
                if ~ismember(FUS_label(r, c), ind_FUS)
                    FUS_filter(r, c) = 0;
                end
            end
        end
        
        %% check for good granule masking
        maskRGB = cat(3, DAPI_filter, clean_granules, FUS_filter);
        figure, montage({origRGB, maskRGB})
        
        Y = 1; N = 0; y = 1; n = 0;
        good_masks = input('Masks okay? (Y or N) ');
        while good_masks ~= 1 && good_masks ~= 0
            good_masks = input('Masks okay? (Y or N) ');
        end
        
        %% if first mask bad, try alternative granule mask
        if ~good_masks
            norm = uint8(255*mat2gray(FUS_o));
            adjust = imadjust(norm, [0.2 1], [0 1], 1.3);
            thresh = round(mean(adjust,'all'));
            BG = (adjust > thresh).*thresh;
            blur = imgaussfilt(norm, 5);
            
            filter = adjust - blur - uint8(BG);
            filter(filter < 0) = 0;
            BW = imbinarize(filter);
            BW2 = bwmorph(BW, 'majority', 3);
            %clean_granules = bwareaopen(BW2, 10);
            clean_granules = BW2;
            
            maskRGB = cat(3, DAPI_filter, clean_granules, FUS_filter);
            figure, montage({origRGB, maskRGB})
            
            good_masks = input('Alternative masks okay? (Y or N) ');
            while good_masks ~= 1 && good_masks ~= 0
                good_masks = input('Invalid input. Alternative masks okay? (Y or N) ');
            end
            
            % update granule mask method used
            if good_masks
                gran_mask_method = 4;
            end
        end
        
        %% analysis
        if good_masks
            %% remove granules in nuclei and not in FUS_filter
            % granule and nucleus overlap
            gran_nuc_overlap = clean_granules + (DAPI_filter > 0) > 1;
            
            gran_label = bwlabel(clean_granules);
            gran_in_FUS = unique(gran_label(FUS_filter >= 1));
            
            if remove_nuc_gran % remove granules with 100% overlap with nuclei
                gran_in_cyto = unique(gran_label(gran_nuc_overlap == 0));
                ind_gran = intersect(gran_in_FUS, gran_in_cyto);
            else
                ind_gran = gran_in_FUS;
            end
            
            if ind_gran(1) == 0
                ind_gran = ind_gran(2:end); % remove zero index
            end
            
            gran_filter = clean_granules;
            for r = 1:size(gran_label, 1)
                for c = 1:size(gran_label, 2)
                    if ~ismember(gran_label(r, c), ind_gran)
                        gran_filter(r, c) = 0;
                    end
                end
            end
            
            %% final masks
            finalRGB = cat(3, DAPI_filter, gran_filter, FUS_filter);
            %figure, imshow(finalRGB)
            %title('Final Masks')
            
            %figure, montage({origRGB, finalRGB})
            
            
            %% determine distance from nuclei
            % update labeled masks
            [DAPI_filt_label, num_nuclei] = bwlabel(DAPI_filter);
            [gran_filt_label, num_gran] = bwlabel(gran_filter);
            
            % determine center of granules (regionprops)
            gran_props = regionprops('table', gran_filt_label, DAPI_o, 'Area',...
                'Centroid', 'MajorAxisLength', 'MinorAxisLength',...
                'WeightedCentroid', 'MaxIntensity',...
                'MeanIntensity', 'MinIntensity');
            
            if ~isempty(gran_props) % check if granules detected
                gran_props{:,1} = double(gran_props{:,1});
                
                gran_centers = round(gran_props{:,'WeightedCentroid'}); % use weighted centroid
                gran_x = gran_centers(:,1);
                gran_y = gran_centers(:,2);
                
                %num_gran = length(gran_x);
                
                % figure, imshow(gran_filter)
                % hold on, plot(gran_x, gran_y, 'm*')
                
                % determine coords of nuclei boundaries
                nuc_border = bwperim(DAPI_filt_label > 0); % mask to boundary
                [mesh_x, mesh_y] = meshgrid(1:size(nuc_border,1), 1:size(nuc_border,2));
                nuc_coords = [mesh_x(:), mesh_y(:), nuc_border(:)]; % coords of all pixels
                nuc_coords = nuc_coords(nuc_coords(:,3) == 1,:); % coords of nuc border pixels
                
                nuc_x = nuc_coords(:,1);
                nuc_y = nuc_coords(:,2);
                
                % hold on, plot(nuc_x, nuc_y, 'r.')
                
                % find closest nucleus coord index
                nearest_nuc_ind = zeros(1, num_gran);
                nearest_nuc_dist = zeros(1, num_gran);
                for i = 1:num_gran
                    [min_dist, indMinDist] = distFromBoundary(gran_x(i), gran_y(i), nuc_x, nuc_y);
                    nearest_nuc_ind(i) = indMinDist;
                    nearest_nuc_dist(i) = min_dist;
                end
                
                nearest_nuc = [nuc_x(nearest_nuc_ind), nuc_y(nearest_nuc_ind)];
                % hold on, plot(nearest_nuc(:,1), nearest_nuc(:,2), 'y.')
                
                % use labeled boundaries to isolate closest nucleus
                %figure, imshow(DAPI_filt_label)
                
                nuc_labels = zeros(1, num_gran);
                gran_info = cell(num_gran, 4); % pair granule with mask
                gran_list = []; % to output .txt
                gran_listNF = []; % to output .txt
                check = zeros(512,512);
                
                figure, imshow(finalRGB)
                for i = 1:num_gran
                    gran_listNF = [gran_listNF; gran_props{i, 'Area'}.*PixelSize^2, ...
                        nearest_nuc_dist(i).*PixelSize, DAPI_filt_label(nearest_nuc(i,2), nearest_nuc(i,1))];
                    
                    if nearest_nuc_dist(i) <= max_gran_dist % granules distance filter
                        if gran_props{i, 'Area'} <= max_gran_area % granule area filter
                            nuc_labels(i) = DAPI_filt_label(nearest_nuc(i,2), nearest_nuc(i,1));
                            hold on, plot(nearest_nuc(i,1), nearest_nuc(i,2), 'r*')
                            hold on, plot(gran_x(i), gran_y(i), 'm.')
                            
                            gran_info{i, 1} = [gran_x(i), gran_y(i)].*PixelSize; % center coords
                            gran_info{i, 2} = gran_props{i, 'Area'}.*PixelSize^2; % granule area
                            gran_info{i, 3} = gran_props{i, 'MeanIntensity'}; % mean granule intensity
                            gran_info{i, 4} = DAPI_filt_label == nuc_labels(i); % nearest nuc mask
                            gran_info{i, 5} = gran_filt_label == i; % granule mask
                            gran_info{i, 6} = nuc_labels(i); % nucleus index
                            gran_info{i, 7} = nearest_nuc_dist(i).*PixelSize; % nearest nucleus distance
                            
                            gran_list = [gran_list; gran_info{i,2}, gran_info{i,7}, gran_info{i,6}];
                        end
                    end
                    
                    %check = check + nuc_gran{i,4} + nuc_gran{i,5};
                end
                
                title('Granule Nuclei Pairing');
                %legend('nearest nucleus', 'granule center', 'Location', 'NorthWest')
                %TODO: put text labels at center of each nuclei with its index number AND/OR number of granules for that nucleus
                saveas(gcf, [file(1:3), '_masks', '.png']);
                
                %figure, imshow(check>0)
                
                % number of granules per nuclei
                gran_per_nuc = zeros(1, num_nuclei);
                for i = 1:size(gran_list,1)
                    gran_per_nuc(gran_info{i, 6}) = gran_per_nuc(gran_info{i, 6}) + 1;
                end
                
                num_gran = size(gran_list,1); % update number granules after filter
                
                disp(['Granules: ' num2str(num_gran)]);
                disp(['Nuclei: ' num2str(num_nuclei)]);
                disp('Granules per nucleus:')
                disp(gran_per_nuc)
                disp(['Average granules per nucleus: ', num2str(num_gran/num_nuclei)])
                if ~isempty(gran_list)
                    disp(['Average granule area: ', num2str(mean(gran_list(:,1)))])
                    disp(['Total granule area: ', num2str(sum(gran_list(:,1)))])
                end
                
            else % if no granules detected
                num_gran = 0;
                gran_info = {};
                gran_list = [];
                gran_listNF = [];
                
                disp(['Granules: ' num2str(num_gran)]);
                disp(['Nuclei: ' num2str(num_nuclei)]);
                disp('Granules per nucleus: 0')
                disp(' ')
            end
            
            %% save data, masks, and input parameters to .mat
            save([file(1:3) '.mat'], 'gran_info', 'gran_list', 'DAPI_mask', 'FUS_mask', 'clean_granules', ...
                'min_area_um', 'min_overlap_area_um', 'max_gran_dist_um', 'max_gran_area_um', ...
                'watershed_slicing', 'remove_nuc_gran', 'gran_listNF', 'maskRGB', 'finalRGB')
            
            %% output data to .txt files
            if ~isempty(gran_listNF) % check if granules detected
                % granules list (area, nearest nuc dist and index) with no filters
                gran_tableNF = table(gran_listNF(:,1), gran_listNF(:,2), gran_listNF(:,3), ...
                    'VariableNames', {'Area', 'NearestNucDist', 'NearestNucIndex'});
                writetable(gran_tableNF, [file(1:3), '_gran_list_NoFilters'])
                
                if ~isempty(gran_list) % check if granules remain after filters
                    % granules list with granule area and distance filters
                    gran_table = table(gran_list(:,1), gran_list(:,2), gran_list(:,3), ...
                        'VariableNames', {'Area', 'NearestNucDist', 'NearestNucIndex'});
                    writetable(gran_table, [file(1:3), '_gran_list'])
                    
                    % batch stats
                    fprintf(fileID,'%s,%d,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f\n', file, num_gran, sum(gran_list(:,1)), ...
                        mean(gran_list(:,1)), std(gran_list(:,1)), std(gran_list(:,1))/sqrt(size(gran_list,1)), ...
                        num_nuclei, size(gran_listNF,1), sum(gran_listNF(:,1)), mean(gran_listNF(:,1)), ...
                        std(gran_listNF(:,1)), std(gran_listNF(:,1))/sqrt(size(gran_listNF,1)) );
                else
                    writetable(table([]), [file(1:3), '_gran_list'])
                    
                    fprintf(fileID,'%s,%d,%f,%f,%f,%d,%d,%f,%f,%f\n', file, ...
                        0, 0, 0, 0, num_nuclei, 0, 0, 0, 0);
                end
            else
                writetable(table([]), [file(1:3), '_gran_list'])
                writetable(table([]), [file(1:3), '_gran_list_NoFilters'])
                
                fprintf(fileID,'%s,%d,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f\n', file, ...
                        0, 0, 0, 0, 0, num_nuclei, 0, 0, 0, 0, 0);
            end
        else
            continue
            % TODO: manual override option? (aka user can say if there's 0 granules for badly masked images)
        end
    end
    fileID = fclose(fileID);
catch ME % catch errors so .txt file is always closed
    if fileID ~= 0
        fclose(fileID);
    end
    disp(['Error encountered at ' num2str(ME.stack.line) ': ' ME.identifier])
    disp(ME.message)
    disp('Closed output .txt file')
end