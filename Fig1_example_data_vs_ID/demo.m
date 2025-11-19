%% This program corresponds to the simulation code for Figure 1
% Required data files:'demov5.mat' and 'train_data.mat'

% Author: Ziyuan Luo, Wenxuan Xu
% email:  ziyuan.luo@bit.edu.cn
%         xuwenxuan1209@163.com
% Last revision: Oct-28-2025

%------------- BEGIN CODE --------------
clear
% load data
load('demo.mat');
n=10;
m=3;
p=2;
r=1;
N = size(x_real,2);
len_plot = 1:N;

%% Simulation
Eid1=[];Eid2=[];Efo=[];
N=100;
len_plot = 1:N;
for k=1:100
    % initialization
    x0_id = zeros(n,1);
    x0_id1 = zeros(size(A_id1,1),1);
    x0_id2 = zeros(size(A_id2,1),1);

    x0 = 5*rand(n,1)-10;
    x_real = x0;
    x_id1 = x0_id1;
    x_id2 = x0_id2;

    z0 = zeros(size(Sigma,1),1);
    z_all = z0;

    U_test = 3 + randn(m,N); % simulation input

    % loop
    for i = 1:N-1
        x1 = A*x0 + B*U_test(:,i); % original system
        y0 = C*x0;                 % output of the original system
        x1_id1 = A_id1*x0_id1 + B_id1*U_test(:,i) + K1*(y0 - C_id1 * x0_id1); % 9-node observer identification system
        x1_id2 = A_id2*x0_id2 + B_id2*U_test(:,i) + K2*(y0 - C_id2 * x0_id2); % 10-node observer identification system
        z1 = Sigma * [U_test(:,i);y0;C*x1;z0]; % data-driven functional observer

        x_real = [x_real x1];
        x_id1 = [x_id1 x1_id1];
        x_id2 = [x_id2 x1_id2];
        z_all = [z_all z1];

        x0 = x1;
        x0_id1 = x1_id1;
        x0_id2 = x1_id2;
        z0 = z1;
    end

    y_real = C*x_real;
    z_real = L*x_real;
    y_id1 = C_id1*x_id1;
    z_id1 = L_id1*x_id1;
    y_id2 = C_id2*x_id2;
    z_id2 = L_id2*x_id2;

    Nsum=50;
    sumid1=sum((z_id1(:,Nsum:end)-z_real(:,Nsum:end)).^2);
    sumid2=sum((z_id2(:,Nsum:end)-z_real(:,Nsum:end)).^2);
    sumfo=sum((z_all(1,Nsum:end)-z_real(1,Nsum:end)).^2);
    Sum=sum(z_real(:,Nsum:end).^2);

    avid1=sqrt(sumid1)/sqrt(Sum);
    avid2=sqrt(sumid2)/sqrt(Sum);
    avfo=sqrt(sumfo)/sqrt(Sum);

    Eid1=[Eid1;avid1];
    Eid2=[Eid2;avid2];
    Efo=[Efo;avfo];

end

k1=1;
Eid1=k1*Eid1;
Eid2=k1*Eid2;
Efo=k1*Efo;

%% Plot results
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultTextFontSize', 10);
set(0, 'DefaultLegendFontSize', 10);

colors = [
    0.85 0.33 0.1;
    0.2 0.63 0.18;
    0.25 0.41 0.88;
    0.55 0.27 0.07;
    0.47 0.67 0.19;
    0.93 0.69 0.13;
    0.494 0.184 0.556;
    0.466 0.674 0.188;
    0.301 0.745 0.933;
    0.635 0.078 0.184;
];

% Plot: All x_id1 trajectories + z_id1 + z_real
figure('Name','x_id1 Trajectories vs z_id1 and z_real','Position',[100 100 1000 600]);
hold on; grid on; box on;
legendEntries = {};  % Store legend entries

plot(len_plot, z_id1,'--^', 'Color', [0.92 0.69 0.13], 'LineWidth', 1, 'MarkerFaceColor', [0.92 0.69 0.13]);
legendEntries{end+1} = '$\hat{z}$';

plot(len_plot, z_real,'-^', 'Color', [0.0 0.447 0.741], 'LineWidth', 1, 'MarkerFaceColor', [0.0 0.447 0.741]);
legendEntries{end+1} = '$z$';

for i = 1:size(x_id1,1)
    color_idx = mod(i-1, size(colors,1)) + 1;
    plot(len_plot, x_id1(i,:),'--','Color', colors(color_idx,:), 'LineWidth', 1);
    legendEntries{end+1} = sprintf('$\\hat{x}^{(9)}_{%d}$', i);
end

xlabel('Time Steps', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('9 Nodes State', 'FontSize', 12, 'FontName', 'Times New Roman');
legend(legendEntries, 'Location', 'Best', 'Interpreter', 'latex', 'FontSize', 15, 'FontName', 'Times New Roman');
set(gcf, 'Position', [200 200 250 200]);
hold off;

% Plot: All x_id2 trajectories + z_id2 + z_real
figure('Name','x_id2 Trajectories vs z_id2 and z_real','Position',[200 200 1000 600]);
hold on; grid on; box on;
legendEntries = {};

plot(len_plot, z_id1,'--^', 'Color', [0.92 0.69 0.13], 'LineWidth', 1, 'MarkerFaceColor', [0.92 0.69 0.13]);
legendEntries{end+1} = '$\hat{z}$';

plot(len_plot, z_real,'-^', 'Color', [0.0 0.447 0.741], 'LineWidth', 1, 'MarkerFaceColor', [0.0 0.447 0.741]);
legendEntries{end+1} = '$z$';

for i = 1:size(x_id2,1)
    color_idx = mod(i-1, size(colors,1)) + 1;
    plot(len_plot, x_id2(i,:),'--','Color', colors(color_idx,:), 'LineWidth', 1);
    legendEntries{end+1} = sprintf('$\\hat{x}^{(10)}_{%d}$', i);
end

xlabel('Time Steps', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('10 Nodes State', 'FontSize', 12, 'FontName', 'Times New Roman');
legend(legendEntries, 'Location', 'Best', 'Interpreter', 'latex', 'FontSize', 15, 'FontName', 'Times New Roman');
set(gcf, 'Position', [200 200 250 200]);
hold off;

% Plot: All z_all trajectories + z_real
figure('Name','z_all Trajectories vs z_real','Position',[300 300 1000 600]);
hold on; grid on; box on;
legendEntries = {};

plot(len_plot, z_all(1,:) ,'--^','Color', [0.92 0.69 0.13], 'LineWidth', 1, 'MarkerFaceColor', [0.92 0.69 0.13]);
legendEntries{end+1} = '$\hat{z}$';

plot(len_plot, z_real,"-^", 'Color', [0.0 0.447 0.741], 'LineWidth', 1, 'MarkerFaceColor', [0.0 0.447 0.741]);
legendEntries{end+1} = '$z$';

for i = 2:size(z_all,1)
    color_idx = mod(i-1, size(colors,1)) + 1;
    plot(len_plot, z_all(i,:) ,'--','Color', colors(color_idx,:), 'LineWidth', 1.2);
    legendEntries{end+1} = sprintf('$\\hat{z}_{%d}^{''}$', i);
end

xlabel('Time Steps', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('FO state', 'FontSize', 12, 'FontName', 'Times New Roman');
legend(legendEntries, 'Location', 'Best', 'Interpreter', 'latex', 'FontSize', 15, 'FontName', 'Times New Roman');
set(gcf, 'Position', [200 200 250 200]);
hold off;

% Plot: y_train and U_train
figure('Name','Training Data Trajectories','Position',[100 100 1000 600]);
hold on; grid on; box on;
legendEntries = {};

for i = 1:size(C*X_train,1)
    y_train=C*X_train;
    color_idx = mod(i-1, size(colors,1)) + 1;
    plot(1:50, y_train(i,:),'--','Color', colors(color_idx,:), 'LineWidth', 1.2);
    legendEntries{end+1} = sprintf('${y}_{%d}$', i);
end

for i = 1:size(U_train,1)
    color_idx = mod(i-1, size(colors,1)) + 2;
    plot(1:50, U_train(i,:),'--','Color', colors(color_idx,:), 'LineWidth', 1.2);
    legendEntries{end+1} = sprintf('${u}_{%d}$', i);
end

xlabel('Time Steps', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('data', 'FontSize', 12, 'FontName', 'Times New Roman');
legend(legendEntries, 'Location', 'Best', 'Interpreter', 'latex', 'FontSize', 15, 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
set(gcf, 'Position', [200 200 250 200]);
hold off;

%%
% error
figure(p+2)
plot(len_plot,z_id1(1,len_plot)-z_real(1,len_plot),'--d','Color','blue','MarkerFaceColor',[0,0,1],'LineWidth', 1,'MarkerSize', 2)
hold on
plot(len_plot,z_id2(1,len_plot)-z_real(1,len_plot),'--s','Color','red','MarkerFaceColor',[1,0,0],'LineWidth', 1,'MarkerSize', 2)
hold on
plot(len_plot,z_all(1,len_plot)-z_real(1,len_plot),'--hexagram','Color','green','MarkerFaceColor',[0,1,0],'LineWidth', 1,'MarkerSize', 2)
legend('9-node ID', '10-node ID', 'data-driven');
ylabel('$(\hat z - z ) \rm (error)$','Interpreter','latex','FontWeight', 'bold', 'FontName', 'Arial');
set(gcf, 'Position', [200 200 250 200]); 

set(gca, 'LooseInset', get(gca,'TightInset'));

%%
% Plot: observer order bar chart
systems = {'9-node ID', '10-node ID', 'Data-driven'};
orders = [9, 10, 4];

figure('Color', 'white');

color1 = [222/255, 88/255, 43/255];
color2 = [243/255, 163/255, 50/255];
color3 = [24/255, 104/255, 178/255];

colors = [
    222/255, 88/255, 43/255;
    24/255, 104/255, 178/255;
    1/255, 138/255, 103/255;
    243/255, 163/255, 50/255
];

barWidth = 0.4;
b1 = bar(1, orders(1), barWidth, 'FaceColor', color1, 'EdgeColor', 'black', 'LineWidth', 0.8);
hold on;
b2 = bar(2, orders(2), barWidth, 'FaceColor', color2, 'EdgeColor', 'black', 'LineWidth', 0.8);
b3 = bar(3, orders(3), barWidth, 'FaceColor', color3, 'EdgeColor', 'black', 'LineWidth', 0.8);
hold off;

xlabel('methods', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Observer order', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'XTick', 1:3, 'XTickLabel', systems, ...
    'FontName', 'Arial', 'FontSize', 9, ...
    'TickDir', 'out', 'Box', 'off');

set(gca, 'YLim', [0, max(orders) + 2], ...
    'YTick', 0:2:12, ...
    'FontName', 'Arial', 'FontSize', 9, ...
    'TickDir', 'out');

for i = 1:length(orders)
    text(i, orders(i) + 0.3, num2str(orders(i)), ...
         'HorizontalAlignment', 'center', ...
         'FontName', 'Arial', 'FontSize', 9, ...
         'FontWeight', 'bold');
end

grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', [0.85, 0.85, 0.85], 'GridAlpha', 0.7);
set(gca, 'XGrid', 'off');

axis tight;
set(gcf, 'Position', [100, 100, 600, 400]);
box off;

%% Plot boxplot
data = [(Eid1), (Eid2), (Efo)];

set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultAxesFontSize', 14);
set(groot, 'DefaultTextFontSize', 14);

figure('Position', [100, 100, 600, 450], 'Color', 'w');

boxColors = [0.2 0.4 0.8;
             0.9 0.5 0.1;
             0.2 0.6 0.3];

h = boxplot(data, ...
    'Labels', {'9-node ID', '10-node ID', 'Data-driven'}, ...
    'Notch', 'off', ...
    'MedianStyle', 'line', ...
    'Symbol', 'o', ...
    'OutlierSize', 4, ...
    'Whisker', 1.5, ...
    'Colors', 'k');

boxes = findobj(gca, 'Tag', 'Box');

for i = 1:length(boxes)
    patch(get(boxes(i), 'XData'), get(boxes(i), 'YData'), boxColors(length(boxes)-i+1,:), ...
          'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 1.5);
end

medians = findobj(gca, 'Tag', 'Median');
set(medians, 'Color', 'k', 'LineWidth', 2);

set(findobj(gca,'Tag','Upper Whisker'), 'Color', 'k', 'LineWidth', 1.2);
set(findobj(gca,'Tag','Lower Whisker'), 'Color', 'k', 'LineWidth', 1.2);
set(findobj(gca,'Tag','Upper Adjacent Value'), 'Color', 'k', 'LineWidth', 1.2);
set(findobj(gca,'Tag','Lower Adjacent Value'), 'Color', 'k', 'LineWidth', 1.2);

outliers = findobj(gca,'Tag','Outliers');
set(outliers, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.7 0.7 0.7]);

grid on;

box on;
set(gca, 'LineWidth', 1.2, 'TickDir', 'out', 'TickLength', [0.015 0.015]);

xlabel('methods', 'FontWeight', 'bold');
ylabel('log_{10}(RRMSE)', 'FontWeight', 'bold');

set(gca, 'LooseInset', get(gca,'TightInset'));
set(gca, 'YScale', 'log');
ylabel('RRMSE', 'FontWeight', 'bold');

