%open the folders containing all the csv files that need to display. Rename
%them as "#.xxxxxx.csv". And then run this code.
clear all;
Dir=pwd;
Target=append(Dir,'\','*.csv');
file_list=dir(Target);
no_file=length(file_list);
file_name=string(ones(no_file,1));
pure_names=string(ones(no_file.*2,1));
figure
hold on
for j=1:no_file
    file_name(j)=convertCharsToStrings(file_list(j).name);
    split_name=split(file_name(j),'.');
    pure_name=string(split_name(2));
    pure_names(2*j-1)=pure_name;
    pure_names(2*j)='Fitting';
    A=readmatrix(file_name(j));
    B=readcell(file_name(j));
    Diameter=A(:,1);
    Value=A(:,2);
    PeakDiameters=Diameter(Value~=0);
    PeakValues=Value(Value~=0);
    LogPeakDiameters=log10(PeakDiameters);
    
    Data_to_Fit=[];
    for i=1:length(LogPeakDiameters)
        number=floor(PeakValues(i).*50);
        Add=ones(number,1).*LogPeakDiameters(i);
        Data_to_Fit=[Data_to_Fit;Add];
    end
    
    hBar=bar(LogPeakDiameters, PeakValues,0.6,'FaceAlpha',0.6, 'EdgeColor', 'None');
    
    FitModel=fitgmdist(Data_to_Fit,2,'SharedCov',true);
    disp(pure_name)
    fit_x=[0:0.01:5]';
    fit_y=pdf(FitModel,fit_x);
    fit_y=fit_y.*10;
    
    plot(fit_x, fit_y, 'linewidth',2, 'color', get(hBar,'facecolor'));
end
%set(gcf, 'Color', 'None')
%set(gca, 'Color', 'None')
set(gca,'Xtick',1:5);
set(gca,'Ylim',[0 100]);
set(gca,'Xticklabel',10.^get(gca,'Xtick'));
xlabel('Diameter/nm');
ylabel('Mass%');
legend(pure_names);
hold off
% filename='E8-6min.csv';
% A=readmatrix(filename);
% B=readcell(filename);
% Diameter=A(:,1);
% Value=A(:,2);
% PeakDiameters=Diameter(Value~=0);
% PeakValues=Value(Value~=0);
% LogPeakDiameters=log(PeakDiameters);
% hBar=bar(LogPeakDiameters, PeakValues,1);
% set(gca,'Xtick',-1:5);
% set(gca,'Xticklabel',10.^get(gca,'Xtick'));

