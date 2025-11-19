colors = [136,86,167;
    49,130,189]/255;
load('./error_observable_ER.mat',"Y1","Y2");
data2 = [log10(Y1);log10(Y2)];
group_inx = [ones(1,size(Y1,1)), 2.*ones(1,size(Y2,1))];
load('./error_unobservable_ER.mat',"Y1","Y2");
data2 = [data2;log10(Y1);log10(Y2)];
group_inx = [group_inx,ones(1,size(Y1,1)), 2.*ones(1,size(Y2,1))];
condition_names = {'100', '120', '140', '160', '180', '200'};
h = daviolinplot(data2,'groups',group_inx,'outsymbol','k+',...
    'xtlabels', condition_names,'color',colors,'scatter',2, 'scattersize',5, ...
    'jitter',1,...
    'box',0,'scattercolors','same',...
    'boxspacing',1.2);

xl = xlim; xlim([xl(1)-0.1, xl(2)+0.1]); % make more space for the legend
set(gca,'FontSize',10);
ylim([-10,0.7]);

set(gca,'ytick',-9:4:0);
set(gcf, 'Position', [400 400 640 280]); 
ax = gca;
ax.YColor = [0, 0, 0]; 
ax.XColor = [0, 0, 0]; 
ax.LineWidth = 0.6; 
box on