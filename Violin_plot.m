clear
clear all 
close all 

[file,path] = uigetfile('*.*');
directory=[path file];
[num,txt,raw] = xlsread(directory);
size = size(raw);
xlabel={raw{1,:}};
for i=1:size(1,2)
    Y{:,i}=num(:,i);
end

[h,L,MX,MED,bw,F,raw_data]=violin_mod(Y,'facecolor','w',...
    'edgecolor','k',...
    'bw', 0.025 ,...
    'mc',[],'medc','r');

ylabel('Dynamic Fraction (f)','FontSize',14)
% print(gcf,'V-Plot.png','-dpng','-r300');  