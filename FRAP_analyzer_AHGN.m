% Written by Amirhossein Ghanbari Niaki [August 2017]
% Description for the code is at the end of the code
% Make sure you have the TIFFStack package installed and added to your
%   paths.

close all 
clc 

%% Time regime & Parameters (make sure you adhust these values)
% pre-bleaching steps all in seconds
interval_1=1; dur_1=1; loop_1=2;
interval_2=0; dur_2=0.1; loop_2=1;
bleach_time=5.7;

%post-bleaching steps all in seconds
interval_3=3; dur_3=120; loop_3=41;
interval_4=10; dur_4=480; loop_4=49;

% creating the time points based on the time regime
time=[-(bleach_time+interval_1) -bleach_time 0];
time=[time linspace(interval_3+time(end),interval_3+time(end)+dur_3,loop_3)];
time=[time linspace(interval_4+time(end),interval_4+time(end)+dur_4,loop_4)];

%other initial parameters
shape='circ';
width=1024;  height=1024;

%% Directory
Directory=input('Directory:: ','s');
	if isempty(Directory)
        Directory=pwd;
    end
cd(Directory);
filename=dir([Directory, '\*.tif']);
number_of_files=length(filename);

% generating necessary variables
t_frames=length(time);
spot_int=NaN(1,t_frames);
ref_int=NaN(number_of_files,t_frames);
bg_int=NaN(number_of_files,t_frames);
int_correction=NaN(1,t_frames);
norm_FRAP=NaN(1,t_frames);
final_report=NaN(t_frames,1+number_of_files+3,5);
header={'Time'};

% this is the main loop in the body of the code going over each .tif file
for file_cnt=1:number_of_files 
    img_info=imfinfo(filename(file_cnt).name);
    frames=length(img_info); 
    
    % read the tiff stack image
    tiff_image=TIFFStack(filename(file_cnt).name);
    
    % generating an RGB image for checking the tranlational drift in the
    % image
    qc_rgb_img=uint8(zeros(width,height,3));
    nb_adj=imadjust(tiff_image(:,:,3));  b_adj=imadjust(tiff_image(:,:,frames));
    nb_img=im2uint8(nb_adj);
    b_img=im2uint8(b_adj);
    qc_rgb_img(:,:,2)=nb_img;  qc_rgb_img(:,:,1)=b_img;
    
    % showing the image
    figure;
    ini_img=gcf;
    figure(ini_img);
    imshow(qc_rgb_img);
    
    % drift correction and check
    quality=input('Is the image overlaping? (Ent or n)  ','s');
    if quality == 'n'
        while true
            imshow(qc_rgb_img);
            input('Zoom into a region and press Enter. ')
            disp('Click on corresponding points (GREEN first then RED)')
            [sx,sy]=ginput(2);
            delta_x=sx(2,1)-sx(1,1);      delta_y=sy(2,1)-sy(1,1);
            x_inc=delta_x/time(1,frames);     y_inc=delta_y/time(1,frames);
            close(ini_img);
            new_image(:,:,1)=tiff_image(:,:,1);     new_image(:,:,2)=tiff_image(:,:,2);        new_image(:,:,3)=tiff_image(:,:,3);
            for i=4:frames
                shift_mat=round([-time(1,i)*y_inc -time(1,i)*x_inc]);
                new_image(:,:,i)=circshift(tiff_image(:,:,i),shift_mat);
            end  
            post_qc_rgb_img=uint8(zeros(width,height,3));
            post_nb_adj=imadjust(new_image(:,:,1));  post_b_adj=imadjust(new_image(:,:,frames));
            post_nb_img=im2uint8(post_nb_adj);
            post_b_img=im2uint8(post_b_adj);
            post_qc_rgb_img(:,:,2)=post_nb_img;  post_qc_rgb_img(:,:,1)=post_b_img;
            % showing the image
            figure;
            ini_img=gcf;
            figure(ini_img);
            imshow(post_qc_rgb_img);
            q=input('Is it fine now? (Ent or n)','s');
            if q == 'n'
                continue;
            else
                close(ini_img);
                break;
            end
        end
    else
        close(ini_img);
        for i=1:frames
            new_image(:,:,i)=tiff_image(:,:,i);
        end
    end
 %% ROI Determination 
    % showing and RGB image for ROI determination
    rgb_img=uint8(zeros(width,height,3));
    nb_adj=imadjust(new_image(:,:,2));  b_adj=imadjust(new_image(:,:,3));
    nb_img=im2uint8(nb_adj);
    b_img=im2uint8(b_adj);
    rgb_img(:,:,2)=nb_img;  rgb_img(:,:,1)=b_img;
    neg_rgb_img=imcomplement(rgb_img);
    
    % showing the image
    figure;
    ini_img=gcf;
    figure(ini_img);
    imshow(neg_rgb_img);
    
    % determining the number or stimulated ROIs
    while true
        n_roi=input('How many ROIs are being analyzed? ');
        if isempty(n_roi)
            continue;
        else
            break;
        end
    end
    % stalling till the zooming is done
    input('ZOOM IN TO THE DESIRED REGION');
    if strcmp(shape,'circ')
        th = 0:pi/50:2*pi;
        s_ROI=zeros(length(th),2*n_roi);
        % finding the bleach spot ROI
        n_cnt=1;
        while n_cnt<=n_roi
            disp('Bleached ROI: Choose the center and the edge')
            [sx,sy]=ginput(2);
            s_center=[sx(1,1) sy(1,1)];
            s_radius=sqrt((sx(2,1)-sx(1,1))^2+(sy(2,1)-sy(1,1))^2);
            xunit = s_radius * cos(th) + sx(1,1);
            yunit = s_radius * sin(th) + sy(1,1);
            hold on
            spot=plot(xunit, yunit,'b','LineWidth',1);
            hold off
            quality_check = input('Is it OK? (Ent or n) ','s');
            if quality_check == 'n'
                set(spot,'Visible','off')
                continue;
            else
                % the coordinates for the spot (x1,y1,x2,y2,...)
                s_ROI(:,2*n_cnt-1:2*n_cnt)=[xunit' yunit'];
                n_cnt=n_cnt+1;
            end
        end

        % finding the Reference ROI
        while true
            disp('Reference ROI: Choose the center and edge')
            [rx,ry]=ginput(2);
            r_center=[rx(1,1) ry(1,1)];
            r_radius=sqrt((rx(2,1)-rx(1,1))^2+(ry(2,1)-ry(1,1))^2);
%                 r_radius=s_radius;
            xunit = r_radius * cos(th) + rx(1,1);
            yunit = r_radius * sin(th) + ry(1,1);
            hold on
            spot=plot(xunit, yunit,'r','LineWidth',1);
            hold off
            quality_check = input('Is it OK? (Ent or n) ','s');
            if quality_check == 'n'
                set(spot,'Visible','off')
                continue;
            else
                r_ROI=[xunit' yunit'];              
            end
            break;
        end

        % finding the Background ROI
        while true
            disp('Backgroung ROI: Choose the center and edge')
            [bx,by]=ginput(2);
            b_center=[bx(1,1) by(1,1)];
            b_radius=sqrt((bx(2,1)-bx(1,1))^2+(by(2,1)-by(1,1))^2);
%                 b_radius=s_radius;
            xunit = b_radius * cos(th) + bx(1,1);
            yunit = b_radius * sin(th) + by(1,1);
            hold on
            spot=plot(xunit, yunit,'k','LineWidth',1);
            hold off
            quality_check = input('Is it OK? (Ent or n) ','s');
            if quality_check == 'n'
                set(spot,'Visible','off')
                continue;
            else
                b_ROI=[xunit' yunit'];              
            end
            break;
        end
        close(ini_img)
    end
       
%% Intensity Measurements
    % making binary masks of ROIs and measuring the mean intensity in the ROI through all frames
    temp_spot_int=NaN(n_roi,t_frames);
    for i=1:frames
        image=new_image(:,:,i);
        for spot_cnt=1:n_roi
            temp_spot_int(spot_cnt,i)=mean(image(poly2mask(s_ROI(:,2*spot_cnt-1),s_ROI(:,2*spot_cnt),width,height)));
        end
        ref_int(file_cnt,i)=mean(image(poly2mask(r_ROI(:,1),r_ROI(:,2),width,height)));
        bg_int(file_cnt,i)=mean(image(poly2mask(b_ROI(:,1),b_ROI(:,2),width,height)));
    end
    spot_int=cat(1,spot_int,temp_spot_int);
    % correcting the intesities based on the background and the reference
    background=mean(bg_int(file_cnt,1:frames));
    temp_int_corr=NaN(n_roi,t_frames);
    temp_norm_frap=NaN(n_roi,t_frames);
    for spot_cnt=1:n_roi
        % correcting for the background and the bleaching
        temp_int_corr(spot_cnt,1:frames)=(temp_spot_int(spot_cnt,1:frames)-background)./(ref_int(file_cnt,1:frames)-background);
        %normalizing based on the pre-bleached frame(2nd frame) and the post-bleached frame(3rd frame)
        temp_norm_frap(spot_cnt,1:frames)=(temp_int_corr(spot_cnt,1:frames)-temp_int_corr(spot_cnt,3))./(temp_int_corr(spot_cnt,2)-temp_int_corr(spot_cnt,3));
        % generating the header
        header{1,end+1}=sprintf('%s_spot%d',filename(file_cnt).name,spot_cnt);
    end
    int_correction=cat(1,int_correction,temp_int_corr);
    norm_FRAP=cat(1,norm_FRAP,temp_norm_frap);
    clear tiff_image
end

%% Post-prcessing and Exporting Data
% getting rid of NaN in the frist line of the data

% measuring the average and std of each time point
data_size=size(spot_int);
norm_FRAP(end+1,:)=nanmean(norm_FRAP);  norm_FRAP(end+1,:)=nanstd(norm_FRAP(1:data_size(1,1),:));  norm_FRAP(end+1,:)=norm_FRAP(end,:)./sqrt(sum(~isnan(norm_FRAP(1:data_size(1,1),:))));
int_correction(end+1,:)=nanmean(int_correction);  int_correction(end+1,:)=nanstd(int_correction(1:data_size(1,1),:));  int_correction(end+1,:)=int_correction(end,:)./sqrt(sum(~isnan(int_correction(1:data_size(1,1),:))));
spot_int(end+1,:)=nanmean(spot_int);  spot_int(end+1,:)=nanstd(spot_int(1:data_size(1,1),:));  spot_int(end+1,:)=spot_int(end,:)./sqrt(sum(~isnan(spot_int(1:data_size(1,1),:))));
ref_int(end+1,:)=nanmean(ref_int);  ref_int(end+1,:)=nanstd(ref_int(1:number_of_files,:));  ref_int(end+1,:)=ref_int(end,:)./sqrt(sum(~isnan(ref_int(1:number_of_files,:))));
bg_int(end+1,:)=nanmean(bg_int);  bg_int(end+1,:)=nanstd(bg_int(1:number_of_files,:));  bg_int(end+1,:)=bg_int(end,:)./sqrt(sum(~isnan(bg_int(1:number_of_files,:))));

% preparing the excel sheet to be exported
ext=input('What to add to the name: ','s');
xls_filename=sprintf('Summary_Report_%s.xlsx',ext);
header{1,end+1}='mean'; header{1,end+1}='std'; header{1,end+1}='se';
xlswrite(xls_filename,header,1);
xlswrite(xls_filename,[time' norm_FRAP(2:end,:)'],1,'A2'); 
xlswrite(xls_filename,[time' int_correction(2:end,:)'],2,'A2');
xlswrite(xls_filename,[time' spot_int(2:end,:)'],3,'A2');
xlswrite(xls_filename,[time' ref_int'],4,'A2');
xlswrite(xls_filename,[time' bg_int'],5,'A2');

% plot the normalized data
figure;
for i=1:file_cnt
    plot(time(1,:),ref_int(i,:),'-');
    hold on
end
hold off