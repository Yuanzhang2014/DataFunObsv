
%: Performance of the observer obtained from the noisy data via Theorem 4
% in Example 1
% History input and output data used in Example 1

                  
% Author: Yuan Zhang, Wenxuan Xu
% email:  xuwenxuan1209@163.com
% Last revision: Oct-28-2025

%------------- BEGIN CODE --------------

clear
close all;
% Estimation with noisy history output data. The truncated SVD is based on
% the singular values of 1/TDD'.
% Example
A=[-2 1 0;1 -3 1;0 0 -1]/4;
C=[1 0 0;0 1 0];
B=[1;2;1]
%L=[1 0 1;0 1 0];
L=[1 0 1];
N=200;
sigma=0.5;
%sigma=0;
times=50;
%sigma is the standard deviation of the entry-wise measurement noise
% data_driven_min_funobserver_measurenoise(a, b, c, L, 10, 0.5, 50)

%times 
%%%%%%%%%%%%%%%%%%%%%%%%% System matrix parameters
% Calculate system dimension information
n = size(A, 1);
m = size(B, 2);
p1 = size(C, 1);
r1 = size(L, 1); % Get the number of rows of L
a = A;
b = B;
c = C;
L0=L;
%  sigma=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%% History data collection
% Generate input data Up with length determined by N, random integers from 1 to 2
Up = randi(2, 1, N);
% Up = randn(1, N);
x0 = rand(n, 1);
% Recursively calculate state variables x1 to x(N-1)
x = x0;
for i = 1:N - 1
    x(:, i + 1) = a * x(:, i) + b * Up(i);
end
% Organize relevant data vectors
data_matrix = [Up; x];
rank(data_matrix);
Up(N + 1) = rand(1);
up = Up(1:N-1);
xp = x(:, 1:end - 1);
yp = c * xp;
yp0=yp;
xf = x(:, 2:end);
yf = c * xf;
yf0=yf;
%% adding measurement noises
yfnoise=normrnd(0,sigma,size(yp));
yp=yp+yfnoise;
yf(:,1:end-1)=yp(:,2:end);
yf(:,end)=yf(:,end)+normrnd(0,sigma,size(yp,1),1);

zp = L * xp;
zf = L * xf;

%%%%%%%%%%%%%%%%%% observer design
p = [up; yp; yf; zp];
f = [p; zf];
% Sigma = zf * (pinv(p));
%%% corrected by noise
T=N-1;
correctI=blkdiag(zeros(size(up,1)),eye(size(yp,1))*sigma^2,eye(size(yf,1))*sigma^2,zeros(size(zp,1)));
p_correct=(p*p'-(N-1)*correctI*sigma^2)/(N-1);
%svd_p= svd(p)/sqrt(T);
svd_p= svd(p*p'/T);
svd_trun_index=find(svd_p>sigma^2*1.1);
[U1,S1,V1] = svd(p);
ptrunc=U1(:,svd_trun_index)*S1(svd_trun_index,svd_trun_index)*V1(:,svd_trun_index)';
Sigma = zf * ptrunc'*pinv(p_correct)/T;
Sigmazp = Sigma(:, (m + p1 + p1 + 1):end);


[U2,~,~] = svd(p_correct);
nu=U2(:,svd_trun_index+1:end)';
if size(nu,1)>0
    %    nu = null(p')';
    nuk = nu(:, (m + p1 + p1 + 1):end);
    desired_poles = 0.4 * rand(size(Sigmazp, 1), 1);
    KK = place(Sigmazp', nuk', desired_poles);
    eig(Sigmazp - KK' * nuk);
    Sigma = Sigma - KK' * nu;
    Sigmazp = Sigma(:, (m + p1 + p1 + 1):end);
    eig(Sigmazp);
end


disp('Order of the functional observer');
disp(size(Sigmazp,1));
d = 40;
X0 = xf(:, end);
U = 0;
Y0 = c * X0; Z0=L*X0;
hZ0 = zeros(r1, 1) - rand(r1, 1)*(0.2);
%     Z0 = zeros(r1, 1);
hZall = [];eZall=[];
hZall = [hZall hZ0];eZall=[eZall hZ0-Z0];
Zall = [];
Zall = [Zall, L * X0];


for t = 1:d
    % U = randi(8);
    U=1;
    X1 = A * X0 + B * U;
    Y1 = C * X1;
    Z1 = L * X1;
    Zall = [Zall Z1];
    hZ1 = Sigma * [U; Y0; Y1; hZ0];
    hZ0 = hZ1;
    X0 = X1;
    Y0 = Y1;
    hZall = [hZall hZ1];
    eZall=[eZall hZ1-Z1];
end
%%%%%%% measurement noise 
ynoise=[yp,yf(:,end)];
ytrue=[yp0,yf0(:,end)];
e_num_rows = size(L0, 1);

%%

figure
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultTextFontSize', 10);
set(0, 'DefaultLegendFontSize', 10);


c_true = [0.71 0.51 1.0];   
c_est  = [0.68 0.85 0.90];  

for i = 1:e_num_rows
    subplot(e_num_rows, 1, i); 
    

    plot(1:size(Zall, 2), hZall(i,:), '-o', ...
        'Color',c_est, 'MarkerSize', 3, 'LineWidth', 1); 
    hold on;
    

    plot(1:size(Zall, 2), Zall(i,:), '-^', ...
        'Color', 'b', 'MarkerSize', 3, 'LineWidth', 1);
    
%     xlabel('{Time}', 'FontSize', 12);
    legend({['$\hat z_{', num2str(i), '}(t)$'], ['$z_{', num2str(i), '}(t)$']}, ...
        'Interpreter', 'latex', 'Location', 'best', 'Box','off');
    

    xlim([0 41]);  
    
    grid off; box on;
    set(gcf, 'Position', [200 200 450 330]); 
end






figure
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultTextFontSize', 10);
set(0, 'DefaultLegendFontSize', 10);
e_num_rows = size(L0, 1);
for i = 1:e_num_rows
    subplot(e_num_rows, 1, i); 
    plot(1:size(eZall, 2), eZall(i, :),'b-', 'LineWidth', 1);
%     xlabel('{Time}', 'Fontsize', 12);
%     ylabel({['$e_{', num2str(i), '}(t)$']}, 'Interpreter', 'latex', 'Fontsize', 12);
    set(gcf, 'Position', [200 200 450 330]);
%     if i == 1
%         title('Errors in a data-driven functional observer', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'b');
%     end
    xlim([0 41]);  
    hold on;
end
figure;
hold on;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultTextFontSize', 10);
set(0, 'DefaultLegendFontSize', 10);
T = size(ytrue,2);


col_red   = [0.85 0.33 0.10]; % 
col_blue  = [0.00 0.45 0.74]; % 
col_green = [0.47 0.67 0.19]; % 

% Up: 
plot(1:T, Up(1:end-1), 'o', 'Color',col_red, ...
    'MarkerFaceColor','none', 'DisplayName','true history input data $u(t)$');

% ytrue(1,:): 
plot(1:T, ytrue(1,:), 'o', 'Color',col_blue, ...
    'MarkerFaceColor','none', 'DisplayName','true history output data $y_{1}(t)$');

% ynoise(1,:): 
plot(1:T, ynoise(1,:), 'o', 'Color',col_blue, ...
    'MarkerFaceColor',col_blue, 'DisplayName','history output data $y_{1}(t)$ with noisy');

% ytrue(2,:): 
plot(1:T, ytrue(2,:), 'o', 'Color',col_green, ...
    'MarkerFaceColor','none', 'DisplayName','true history output data $y_{2}(t)$');

% ynoise(2,:):
plot(1:T, ynoise(2,:), 'o', 'Color',col_green, ...
    'MarkerFaceColor',col_green, 'DisplayName','history output data $y_{2}(t)$ with noisy');

legend('Location','best','Interpreter','latex');

ylim([-5,6]);
grid off;
set(gcf, 'Position', [200 200 500 400]); 
hold off;
