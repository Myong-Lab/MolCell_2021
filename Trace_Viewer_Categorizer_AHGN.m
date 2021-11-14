%%%% Original code by Digvijay Singh Edited by Amirhossein Ghanbari Niaki
%%%% READ ME FIRST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%% Trace viewer let's you analyze single molecule traces while catergorizing
%%%%      the traces in different categories for statistical purposes
%%%% I tried to comment as much as possible to make understanding the code easier 
%%%% NOTE: If you need to stop the code and the data be saved Press e instead
%%%%       of ctrl+C
%%%% Set the parameters first and make sure they match your experiments

close all; 
clear all; 


% enter the directory which needs to be analyzed and contains the .trace
% files
Directory_of_TracesFiles=input('Directory:: ','s');
	if isempty(Directory_of_TracesFiles)
        Directory_of_TracesFiles=pwd;
    end
    
    cd(Directory_of_TracesFiles);

%% Set Parameters 
TimeUnit=0.1;
GammaFactor=1.0;
ChannelLeakage=0.2; 
offset=400; %offset


axis_start=1;

min_Int=150; %minimum total intesity 
Cat_Num = 7; % number of categories
spot_viewer='OFF'; % to have the spots shown on the figure turn in ON
green_only='OFF'; % to show PIFE or single green channel traces index that needs to be analyzed (e.g. enter 20 for hel20)
FileIndexNumber=input('Index:: ');
	if isempty(FileIndexNumber)
        FileIndexNumber=1;
    end
GenericFileType='.traces';    
fprintf('Analyzing hel%d%s\n',FileIndexNumber,GenericFileType);

%% reading the chosen .trace file and extractin trace info
File_id=fopen(['hel' num2str(FileIndexNumber) GenericFileType],'r');
Length_of_the_TimeTraces=fread(File_id,1,'int32');

Number_of_traces=fread(File_id,1,'int16');
fprintf('The number of traces in this file is: %d\n',Number_of_traces/2);
% this reads all the raw data into a single matrix appending info of each
% frame
Raw_Data=fread(File_id,Number_of_traces*Length_of_the_TimeTraces,'int16');
disp('Done reading data');
fclose(File_id);

%% Preparing the data for analysis
Index_of_SelectedSpots=(1:Number_of_traces*Length_of_the_TimeTraces);
DataMatrix=zeros(Number_of_traces,Length_of_the_TimeTraces);
Donors=zeros(Number_of_traces/2,Length_of_the_TimeTraces);
Acceptors=zeros(Number_of_traces/2,Length_of_the_TimeTraces);
DataMatrix(Index_of_SelectedSpots)=Raw_Data(Index_of_SelectedSpots);
% assigning the Donor and Acceptor.  
for i=1:(Number_of_traces/2)
   Donors(i,:)=DataMatrix(i*2-1,:);   %So this will be a matrix where each column will be the Donor time series of each selected spot of the movie
   Acceptors(i,:)=GammaFactor.*DataMatrix(i*2,:); %So this will be a matrix where each column will be the Acceptor time series of each selected spot of the movie
end
TimeSeries=(0:(Length_of_the_TimeTraces-1))*TimeUnit;
% reading the peak positions 

pks_id=fopen(['hel' num2str(FileIndexNumber) '.pks'],'r');
pks_info = fscanf(pks_id,'%f %f %f %f',[4 Inf]);
pks_info = pks_info';
fclose(pks_id);
%% Organizing outputs
New_Folder=sprintf('Selected hel%d Traces',FileIndexNumber);
if exist(New_Folder,'dir')==0
    mkdir(New_Folder);
end
TracesCounter=0;
figure;
hd13=gcf;
DT1=[];DT2=[];DT3=[];
DT1a=[];DT2a=[];DT3a=[];
DT1d=[];DT2d=[];DT3d=[];
DT1f=[];DT2f=[];DT3f=[];

Type_counter_cat=zeros(1,Cat_Num);
fraction_list=[];
DT_cnt=0;

%% Notations

% print the options menu
disp('Press Enter to advance');
disp('Press s to Save');
disp('Press b to go Back');
disp('Press bg to subtract the BackGround');
disp('Press c to Cut out a trace');
disp('Press d to measure Dwell times');
disp('Press f to measure Trace Fractions');
disp('Press e to Exit with saving everything');
disp('Press corr to measure the correlation factor')
disp('Enter a number identifying your category');


%% Body of the code 
while TracesCounter < Number_of_traces/2 
   TracesCounter = TracesCounter + 1 ;
   figure(hd13);
   
   % Ploting the donor and acceptor traces
   ax1=subplot(3,4,[1 2 3 4]);
   if strcmp(green_only,'OFF')
       plot(TimeSeries(axis_start:end),Donors(TracesCounter,axis_start:end),'g',...
           TimeSeries(axis_start:end),(Acceptors(TracesCounter,axis_start:end)-ChannelLeakage*Donors(TracesCounter,axis_start:end)),'r',...
           TimeSeries(axis_start:end),offset+Donors(TracesCounter,axis_start:end)+(Acceptors(TracesCounter,axis_start:end)-ChannelLeakage*Donors(TracesCounter,axis_start:end)),'k');
       Trace_ymax = offset+max(Donors(TracesCounter,axis_start:end))+(max(Acceptors(TracesCounter,axis_start:end))-ChannelLeakage*max(Donors(TracesCounter,axis_start:end)));
       axis([axis_start*TimeUnit TimeSeries(end) -100 Trace_ymax+50])
       
       title(sprintf('Molecule %d/%d of hel%d.traces',TracesCounter,Number_of_traces/2,FileIndexNumber));
       
       % Ploting the FRET trace
       ax2=subplot(3,4,[5 6 7 8]);
       FRET_Time_Series=(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:))...
           ./(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:)+(Donors(TracesCounter,:)));

       for m=1:Length_of_the_TimeTraces % getting rid of timepoints with lower intesities than minInt
           if Acceptors(TracesCounter,m)+Donors(TracesCounter,m)<min_Int
               FRET_Time_Series(m)=NaN;
           end
       end
       plot(TimeSeries(axis_start:end),FRET_Time_Series(axis_start:end),'LineWidth',0.5,'Color','b');
       axis([axis_start*TimeUnit TimeSeries(end) -0.1 1])

       linkaxes([ax1,ax2],'x');
   else
       plot(TimeSeries(axis_start:end),Donors(TracesCounter,axis_start:end),'g')
       Trace_ymax = max(Donors(TracesCounter,axis_start:end));
       axis([axis_start*TimeUnit TimeSeries(end) -500 Trace_ymax+50])
       title(sprintf('Molecule %d/%d of hel%d.traces',TracesCounter,Number_of_traces/2,FileIndexNumber));
       
       
       ax2=subplot(3,4,[5 6 7 8]);
       histogram(Donors(TracesCounter,axis_start:end),25,'Normalization','probability','BinMethod','auto')
       
       FRET_Time_Series=(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:))...
           ./(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:)+(Donors(TracesCounter,:)));

       for m=1:Length_of_the_TimeTraces % getting rid of timepoints with lower intesities than minInt
           if Acceptors(TracesCounter,m)+Donors(TracesCounter,m)<min_Int
               FRET_Time_Series(m)=NaN;
           end
       end
   end

   % showing the spot alongside the traces
   if strcmp(spot_viewer,'ON')
       % ploting the FRET histogram of the trace
       subplot(3,4,[9 10]);
       hist(FRET_Time_Series,-.1:.025:1.1);
       xlabel('FRET'); ylabel('Count'); 
       temp=axis;   temp(1)=-0.1;   temp(2)=1.1;   axis(temp);

       % reading the image of the spots from the tiff image
       [X,map] = imread(['hel' num2str(FileIndexNumber) '_ave.tif']);
       if ~isempty(map)
           Im = ind2rgb(X,map);
       end

       % getting the coordinates of the spots from the .pks file
       dspotx = pks_info(TracesCounter*2-1,2); %xcoord
       dspoty = pks_info(TracesCounter*2-1,3); %ycoord
       aspotx = pks_info(TracesCounter*2,2); %xcoord
       aspoty = pks_info(TracesCounter*2,3); %ycoord
       circleradius = 4.5; %dimension of the indicator circle

       % drawing the circle for the spots
       circles = int32([dspotx dspoty circleradius;aspotx aspoty circleradius]);
       donorInserter = insertShape(Im, 'circle', [dspotx dspoty circleradius], 'Color', 'cyan');
       acceptorInserter = insertShape(Im, 'circle', [aspotx aspoty circleradius], 'Color', 'cyan');

       sz=512;
       subplot(3,4,11);
       imshow(imresize(acceptorInserter(int16(max(1,aspoty-20)):int16(min(sz,aspoty+20)),int16(max(1,aspotx-20)):int16(min(512,aspotx+20)),:),4,'nearest'));
       title('Acceptor');
       subplot(3,4,12);
       imshow(imresize(donorInserter(int16(max(1,dspoty-20)):int16(min(sz,dspoty+20)),int16(max(1,dspotx-20)):int16(min(512,dspotx+20)),:),4,'nearest'));
       title('Donor');
       zoom on;
   else
       % ploting the FRET histogram of the trace
       subplot(3,4,[10 11]);
       hist(FRET_Time_Series,-.1:.025:1.1);
       xlabel('FRET'); ylabel('Count'); 
       temp=axis;   temp(1)=-0.1;   temp(2)=1.1;   axis(temp);
   end
   
   %% Choosing what to do next after showing the traces
   
   choice = input(sprintf('molecule#%d; What to do? ',TracesCounter),'s');

   %% Categorization
   cat_id = str2double(choice);
   if ~isnan(cat_id)
       Type_counter_cat(end+1,cat_id) = TracesCounter;
   end

   %% exits the code with saving the categories instead of ctrl+C
   if choice == 'e'
       break;
   end
   
   %% Correlation factor
   if ismember(choice,{'corr','cor','crr','cr','co'})
       [x,y]=ginput(2);
       corr_start = floor((x(1)-TimeSeries(1))/TimeUnit);
       corr_end = floor((x(2)-TimeSeries(1))/TimeUnit);
       if corr_start < corr_end
           corr_factor = corrcoef(Donors(TracesCounter,corr_start:corr_end),Acceptors(TracesCounter,corr_start:corr_end));
       else 
           corr_factor = corrcoef(Donors(TracesCounter,corr_end:corr_start),Acceptors(TracesCounter,corr_end:corr_start));
       end
       resp = sprintf('The Correlation Coefficient is %g', corr_factor(1,2));
       disp(resp)
       TracesCounter=TracesCounter - 1;
       continue;
   end
   
   %% Fraction measurments
   %this option allows the user to choose portions of the trace and
   % calculate the fraction of a behavior in a large fraction of the trace
   if choice == 'f'
       disp('Click twice to indicate your range for calculation');
       [x,y]=ginput(2);
       total_range=abs(x(2)-x(1));
       fraction_range=0;
       Done_Cutting_Choice='n';
       Cut_Counter=0;
       disp('Click twice to indicate your desired fraction range');       
       while true
            [time,y]=ginput;
            seld=size(time);
            hold on
            
            cross=plot(time, y, 'x', 'Color', 'r');
            if mod(seld(1),2) == 0
                temp_choice = input(sprintf('You selected %d points\nEnter to accept and go to next image\nPress r to re-select\n', seld(1)),'s');
                if temp_choice == 'r'
                    set(cross,'Visible','off');  
                    continue
                elseif temp_choice == ''
                    DT_cnt=DT_cnt+(seld(1)/2);
                    continue
                end
            else
                disp('NOT ENOUGH POINTS!!!')
                temp_choice = input('Press r to re-select\nPress e to Exit\n','s');
                if temp_choice == 'r'
                    set(cross,'Visible','off');  
                    continue
                elseif temp_choice == ''
                    continue
                elseif temp_choice == 'e'
                    hold off
                    break
                end
            end
            hold off
            break
       end
       for i=1:2:seld(1)-1
           fraction_range=fraction_range+abs(time(i+1)-time(i));
       end
       fraction_list=[fraction_list (fraction_range/total_range)];
       length(fraction_list)
    end
   
   %% Backwards 
   % this option make the program to go back and show the previous trace
   if choice == 'b'
      TracesCounter=TracesCounter - 2;
      continue;
   end
   
   %% Save traces
   % save the trace as a .dat file with 3 columns; Time, Donor Int, and
   % Corected Acceptor Intensity
   if choice == 's'     
       fname1=[Directory_of_TracesFiles '\' New_Folder '\hel' num2str(FileIndexNumber) '_trace' num2str(TracesCounter) '.dat'];
       Output_to_be_saved_inaFile=[TimeSeries()',Donors(TracesCounter,:)',(Acceptors(TracesCounter,:)-ChannelLeakage*Donors(TracesCounter,:))',FRET_Time_Series'];
       save(fname1,'Output_to_be_saved_inaFile','-ascii') ;
       TracesCounter=TracesCounter - 1;
       continue;
   end
   
   %% single click Save traces
   % save the trace as a .dat file with 3 columns; Time, Donor Int, and
   % Corected Acceptor Intensity
   if strcmp(choice ,'ss')
       disp('Click at the end of the trace');
       Done_Cutting_Choice='n';
       Cut_Counter=0;
       while Done_Cutting_Choice == 'n'
           Cut_Counter=Cut_Counter+1;
           [x,y]=ginput(1);  % Make Two clicks to specify a region which you want to cut out.
           x(1)=round(x(1)/TimeUnit); % ending point of the cut region
           
           
           fname1=[Directory_of_TracesFiles '\' New_Folder '\hel' num2str(FileIndexNumber) '_trace' num2str(TracesCounter) '_truncated.dat'];
           Output_to_be_saved_inaFile=[TimeSeries(1:x(1))' Donors(TracesCounter,1:x(1))' (Acceptors(TracesCounter,1:x(1))-ChannelLeakage*Donors(TracesCounter,1:x(1)))' FRET_Time_Series(1,1:x(1))'];
           save(fname1,'Output_to_be_saved_inaFile','-ascii') ;
           % Asking you whether you want to keep cutting the trace or want
           % to move on to the NEXT set of traces.
           Done_Cutting_Choice=input('Are you done cutting? ("ENT" to finish "n" to continue) ','s');   
           TracesCounter=TracesCounter - 1;
           if Done_Cutting_Choice == 'e'
               break;
           end
           continue;
       end
    end
   
   %% background subtraction
   if strcmp(choice,'bg')
       [x,y]=ginput(2);
       bg_start = floor((x(1)-TimeSeries(1))/TimeUnit);
       bg_end = floor((x(2)-TimeSeries(1))/TimeUnit);
       Cy3_bg=mean(Donors(TracesCounter,bg_start:bg_end));
       Cy5_bg=mean(Acceptors(TracesCounter,bg_start:bg_end));

       Donors(TracesCounter,:) = Donors(TracesCounter,:)-Cy3_bg;
       Acceptors(TracesCounter,:) = Acceptors(TracesCounter,:)-Cy5_bg;  
       TracesCounter=TracesCounter - 1;
       continue;
   end
   
   %% Cut traces 
   % cut a postion of the trace and save it
    if choice == 'c'
       disp('Click twice to select the range which you want to cut out');
       Done_Cutting_Choice='n';
       Cut_Counter=0;
       while Done_Cutting_Choice == 'n'
           Cut_Counter=Cut_Counter+1;
           [x,y]=ginput(2);  % Make Two clicks to specify a region which you want to cut out.
           x(1)=round(x(1)/TimeUnit); % Starting point of the cut region
           x(2)=round(x(2)/TimeUnit); % End point of the cut region
           
           fname1=[Directory_of_TracesFiles '\' New_Folder '\hel' num2str(FileIndexNumber) '_trace' num2str(TracesCounter) '_truncated.dat'];
           Output_to_be_saved_inaFile=[TimeSeries(x(1):x(2))' Donors(TracesCounter,x(1):x(2))' (Acceptors(TracesCounter,x(1):x(2))-ChannelLeakage*Donors(TracesCounter,x(1):x(2)))' FRET_Time_Series(1,x(1):x(2))'];
           save(fname1,'Output_to_be_saved_inaFile','-ascii') ;
           % Asking you whether you want to keep cutting the trace or want
           % to move on to the NEXT set of traces.
           Done_Cutting_Choice=input('Are you done cutting? ("ENT" to finish "n" to continue) ','s');   
           TracesCounter=TracesCounter - 1;
           if Done_Cutting_Choice == 'e'
               break;
           end
           continue;
       end
    end

    %% Dwell time analysis
    if choice == 'd'
        while true
            disp('      Click for beginning and end of states.');
            disp('      Left/middle/right click for different states.');
            [time,y,button]=ginput;
            seld=size(time);
            hold on
            cross=plot(time, y, 'x', 'Color', 'b');
            if mod(seld(1),2) == 0
                temp_choice = input(sprintf('You selected %d points\nEnter to accept and go to next image\nPress r to re-select\n', seld(1)),'s');
                if temp_choice == 'r'
                    set(cross,'Visible','off');  
                    continue
                elseif temp_choice == ''
                    DT_cnt=DT_cnt+(seld(1)/2);
                    continue
                end
            else
                disp('NOT ENOUGH POINTS!!!')
                temp_choice = input('Press r to re-select\nPress e to Exit\n','s');
                if temp_choice == 'r'
                    set(cross,'Visible','off');  
                    continue
                elseif temp_choice == ''
                    continue
                elseif temp_choice == 'e'
                    hold off
                    break
                end
            end
            hold off
            time1=time(button==1);
            for c=1:2:sum(button==1)-1
                t1=ceil(time1(c)/TimeUnit);t2=ceil(time1(c+1)/TimeUnit);
                DT1(end+1)=abs(time1(c+1)-time1(c));
                DT1a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
                DT1d(end+1)=mean(Donors(TracesCounter,t1:t2));
                DT1f(end+1)=mean(FRET_Time_Series(t1:t2));
            end
            time2=time(button==2);
            for c=1:2:sum(button==2)-1
                t1=ceil(time2(c)/TimeUnit);t2=ceil(time2(c+1)/TimeUnit);
                DT2(end+1)=abs(time2(c+1)-time2(c));
                DT2a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
                DT2d(end+1)=mean(Donors(TracesCounter,t1:t2));
                DT2f(end+1)=mean(FRET_Time_Series(t1:t2));
            end
            time3=time(button==3);
            for c=1:2:sum(button==3)-1
                t1=ceil(time3(c)/TimeUnit);t2=ceil(time3(c+1)/TimeUnit);
                DT3(end+1)=abs(time3(c+1)-time3(c));
                DT3a(end+1)=mean(Acceptors(TracesCounter,t1:t2));
                DT3d(end+1)=mean(Donors(TracesCounter,t1:t2));
                DT3f(end+1)=mean(FRET_Time_Series(t1:t2));
            end
            break
        end
    end
    
    %% special dweltime
    if choice == 'q'
        while true
            disp('      Click for the desired location');
            disp('      Left/middle/right click for different states.');
            [time,y,button]=ginput;
            seld=size(time);
            hold on
            cross=plot(time, y, 'x', 'Color', 'b');
            if mod(seld(1),2) == 0
                temp_choice = input(sprintf('You selected %d points\nEnter to accept and go to next image\nPress r to re-select\n', seld(1)),'s');
                if temp_choice == 'r'
                    set(cross,'Visible','off');  
                    continue
                elseif temp_choice == ''
                    DT_cnt=DT_cnt+(seld(1)/2);
                    continue
                end
            else
                disp('NOT ENOUGH POINTS!!!')
                temp_choice = input('Press r to re-select\nPress e to Exit\n','s');
                if temp_choice == 'r'
                    set(cross,'Visible','off');  
                    continue
                elseif temp_choice == ''
                    continue
                elseif temp_choice == 'e'
                    hold off
                    break
                end
            end
            hold off
            time1=time(button==1);
            for c=1:2:sum(button==1)-1
                t1=ceil(time1(c)/TimeUnit);
                DT1(end+1)=time1(c);
                DT1a(end+1)=time1(c+1);
                DT1d(end+1)=abs(time1(c+1)-time1(c));
                DT1f(end+1)=FRET_Time_Series(t1);
            end
            time2=time(button==2);
            for c=1:sum(button==2)
                t1=ceil(time2(c)/TimeUnit);
                DT2(end+1)=time2(c);
                DT2a(end+1)=Acceptors(TracesCounter,t1);
                DT2d(end+1)=Donors(TracesCounter,t1);
                DT2f(end+1)=FRET_Time_Series(t1);
            end
            time3=time(button==3);
            for c=1:2:sum(button==3)-1
                t1=ceil(time3(c)/TimeUnit);
                DT3(end+1)=time3(c);
                DT3a(end+1)=time3(c+1);
                DT3d(end+1)=abs(time3(c+1)-time3(c));
                DT3f(end+1)=FRET_Time_Series(t1);
            end
            break
        end
    end
    
    %% jumping to a trace
    if ismember(choice,{'jump','ju','jp','jum','jmp','jup'})
        Trace_choice=input('-Which # trace do you want to go to?  ');
        TracesCounter=Trace_choice-1;
        continue;
    end
end

%% Post processind the dwell time data
if ~isempty(DT1)
    DT1=[DT1;DT1a;DT1d;DT1f]';
    fname=[Directory_of_TracesFiles '\Dwelltime1.xlsx'];
    if exist('Dwelltime1.xlsx','file') == 2
        temp=xlsread('Dwelltime1.xlsx');
        temp=[temp;DT1];
        xlswrite('Dwelltime1.xlsx',temp);
    end
    xlswrite(fname,DT1)
end
if ~isempty(DT2)
    DT2=[DT2;DT2a;DT2d;DT2f]';
    fname=[Directory_of_TracesFiles '\Dwelltime2.xlsx'];
    if exist('Dwelltime2.xlsx','file') == 2
        temp=xlsread('Dwelltime2.xlsx');
        temp=[temp;DT2];
        xlswrite('Dwelltime2.xlsx',temp);
    end
    xlswrite(fname,DT2)
end
if ~isempty(DT3)
    DT3=[DT3;DT3a;DT3d;DT3f]';    
    fname=[Directory_of_TracesFiles '\Dwelltime3.xlsx'];
    if exist('Dwelltime3.xlsx','file') == 2
        temp=xlsread('Dwelltime3.xlsx');
        temp=[temp;DT3];
        xlswrite('Dwelltime3.xlsx',temp);
    end
    xlswrite(fname,DT3)
end

%% Post processing the category data
max=0;
total_nnz=nnz(Type_counter_cat);
for j=1:Cat_Num
    if nnz(Type_counter_cat(:,j))>max
        max=nnz(Type_counter_cat(:,j));
    end
end
CategoryMatrix=zeros(max,Cat_Num);

for j=1:Cat_Num
    nnz_count=nnz(Type_counter_cat(:,j));
    CategoryMatrix(1:nnz_count,j)=nonzeros(Type_counter_cat(:,j));
    percentage = sprintf('Category #%d: %d (%g%%)',j,nnz_count,(nnz_count/total_nnz)*100);
    disp(percentage)
end
if ~isempty(CategoryMatrix)
    fname=[Directory_of_TracesFiles '\' New_Folder '\Category Report.xlsx'];
    xlswrite(fname,CategoryMatrix);
end
if ~isempty(fraction_list)
    fname=[Directory_of_TracesFiles '\' New_Folder '\Dynamic population.xlsx'];
    xlswrite(fname,fraction_list');
end

%% binding event data export



close all;
















