% Sensor fault detection and recovery in water networks (EPANET’s network 3)
% Author: Yuan Zhang, Ranbo Cheng, Ziyuan Luo
% Email: zhangyuan14@bit.edu.cn
%        ranbo_cheng@bit.edu.cn
%        ziyuan.luo@bit.edu.cn
% Last revision: Oct-20-2025

%------------- BEGIN CODE --------------
% Initialization
% Load data from files 'junction_infoMatrix', 'reservoirs_infoMatrix',
% 'tanks_infoMatrix', 'pipe_infoMatrix', 'pumps_infoMatrix'.
%%
clear;clc;close all;
addpath('../utils')
load('.\pipe_infoMatrix.mat');
load('.\pumps_infoMatrix.mat');
load('.\tanks_infoMatrix.mat');
load('.\junction_infoMatrix.mat');
load('.\reservoirs_infoMatrix.mat');
junc = size(junction_infoMatrix,1);
res = size(reservoirs_infoMatrix,1);
tank = size(tanks_infoMatrix,1);
pip = size(pipe_infoMatrix,1);
pu = size(pumps_infoMatrix,1);
n = junc + tank;
m = 2;
p = 4;
r = 2;
A = zeros(n,n);
B = zeros(n,m);          % Input matrix
B(1,1) = 10;
B(7,2) = 10;
I = eye(n);
C = I([1,11,17,85],:);   % Output matrix
L = I([1,10],:);         % Functional state matrix
dt = 0.01;               % Discrete step size

% Construct the water network system matrix
for i = 1 : pip
    pipe_infoMatrix(i,5)= 10 * rand(1) + 10;       % Pipe diameter
    pipe_infoMatrix(i,4)= 1000 * rand(1) +3000;    % Pipe length
end
for i = 1 : junc 
    for j = 1 : junc
        for l = 1 : pip
            if pipe_infoMatrix(l,2) == junction_infoMatrix(i,1) && pipe_infoMatrix(l,3) == junction_infoMatrix(j,1)
                A(i,i) = A(i,i) - discretized(pipe_infoMatrix(l,5),pipe_infoMatrix(l,4),dt);
                A(i,j) = A(i,j) + discretized(pipe_infoMatrix(l,5),pipe_infoMatrix(l,4),dt);
                A(j,j) = A(j,j) - discretized(pipe_infoMatrix(l,5),pipe_infoMatrix(l,4),dt);
                A(j,i) = A(j,i) + discretized(pipe_infoMatrix(l,5),pipe_infoMatrix(l,4),dt);
            end
        end
        for l = 1 : pu
            if pumps_infoMatrix(l,2) == junction_infoMatrix(i,1) && pumps_infoMatrix(l,3) == junction_infoMatrix(j,1)
                lenPump = 100 * rand(1);
                dPump = 10 * rand(1);
                A(i,i) = A(i,i) - discretized(dPump,lenPump,dt);
                A(i,j) = A(i,j) + discretized(dPump,lenPump,dt);
                A(j,j) = A(j,j) - discretized(dPump,lenPump,dt);
                A(j,i) = A(j,i) + discretized(dPump,lenPump,dt);
            end
        end
    end
end
for i = 1 : tank
    for j = 1 : junc
        for l = 1 : pip
            if pipe_infoMatrix(l,2) == tanks_infoMatrix(i,1) && pipe_infoMatrix(l,3) == junction_infoMatrix(j,1)
                A(i+junc,i+junc) = A(i+junc,i+junc) - discretized(pipe_infoMatrix(l,5),pipe_infoMatrix(l,4),dt);
                A(i+junc,j) = A(i+junc,j) + discretized(pipe_infoMatrix(l,5),pipe_infoMatrix(l,4),dt);
                A(j,j) = A(j,j) - discretized(pipe_infoMatrix(l,5),pipe_infoMatrix(l,4),dt);
                A(j,i+junc) = A(j,i+junc) + discretized(pipe_infoMatrix(l,5),pipe_infoMatrix(l,4),dt);
            end
        end
    end
end
for i = 1 : n
    A(i,i) = A(i,i) + 1;
end

% Building training data
N = n * m + 2 * n + 1;
X = zeros(n,N+1);
X(:,1) = 100 * rand(n,1) + 100 * ones(n,1);     % Randomly generate initial state
U = 100 * rand(2,N) + ones(2,N) * 200;
A = A./max(abs(eig(A))*1.1);
B = B * dt;
for k = 1 : N
    X(:,k+1) = A * X(:,k) + B * U(:,k);
end

% Design observers based on different sensors
[Sigma1,exist1,newL1] = Reduced_order_observ_greedy_trajectory(U,X(:,1:N),C(1,:),L,1);
[Sigma2,exist2,newL2] = Reduced_order_observ_greedy_trajectory(U,X(:,1:N),C(2,:),L,1);
[Sigma3,exist3,newL3] = Reduced_order_observ_greedy_trajectory(U,X(:,1:N),C(3,:),L,1);
[Sigma4,exist4,newL4] = Reduced_order_observ_greedy_trajectory(U,X(:,1:N),C(4,:),L,1);

% Initializing test data
Yt_att = zeros(1,N);
Yt = zeros(p,N);
Zt = zeros(r,N);
Ut = zeros(m,N);
Xt = zeros(n,N+1);
Xt(:,1) = 100 * rand(n,1) + 50 * ones(n,1);    % Randomly generate initial state
Attack = 10 * rand(1,N);                       % Sensor attack value
tt = 161:190;                                  % Sensor attack time period

% Building test data
for k = 1 : tt(1,1)-1
    Ut(:,k) = 300 * ones(m,1) - [Xt(10,k);Xt(60,k)];
    Xt(:,k+1) = A * Xt(:,k) + B * Ut(:,k);
    Yt_att(:,k) = C(1,:) * Xt(:,k);
    Yt(:,k) = C * Xt(:,k);
    Zt(:,k) = L * Xt(:,k);
end
for k = tt(1,1) : tt(1,end)
    Ut(:,k) = 300 * ones(m,1) - [Xt(10,k);Xt(60,k)];
    Xt(:,k+1) = A * Xt(:,k) + B * Ut(:,k);
    Yt_att(:,k) = C(1,:) * Xt(:,k) + Attack(:,k);
    Yt(:,k) = C * Xt(:,k);
    Zt(:,k) = L * Xt(:,k);
end
for k = tt(1,end) + 1 : N
    Ut(:,k) = 300 * ones(m,1) - [Xt(10,k);Xt(60,k)];
    Xt(:,k+1) = A * Xt(:,k) + B * Ut(:,k);
    Yt_att(:,k) = C(1,:) * Xt(:,k);
    Yt(:,k) = C * Xt(:,k);
    Zt(:,k) = L * Xt(:,k);
end

% Observing functional states
Zt_observ_att = zeros(size(Sigma1,1),N);
Zt_observ1 = zeros(size(Sigma1,1),N);
Zt_observ2 = zeros(size(Sigma2,1),N);
Zt_observ3 = zeros(size(Sigma3,1),N);
Zt_observ4 = zeros(size(Sigma4,1),N);
for k = 1 : N - 1
    Zt_observ_att(:,k+1) = Sigma1 * [Ut(:,k);Yt_att(:,k);Yt_att(:,k+1);Zt_observ_att(:,k)];
    Zt_observ1(:,k+1) = Sigma1 * [Ut(:,k);Yt(1,k);Yt(1,k+1);Zt_observ1(:,k)];
    Zt_observ2(:,k+1) = Sigma2 * [Ut(:,k);Yt(2,k);Yt(2,k+1);Zt_observ2(:,k)];
    Zt_observ3(:,k+1) = Sigma3 * [Ut(:,k);Yt(3,k);Yt(3,k+1);Zt_observ3(:,k)];
    Zt_observ4(:,k+1) = Sigma4 * [Ut(:,k);Yt(4,k);Yt(4,k+1);Zt_observ4(:,k)];
end

% Drawing Fig.4(b)
close all;
figure(1);
len_plot = 1:100;
plot(len_plot, Zt(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ1(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ2(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ3(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ4(2,len_plot), 'LineWidth', 1.2);
hold on;
legend('Target functional state', 'FO via Sensor 1','FO via Sensor 2','FO via Sensor 3','FO via Sensor 4')
xlabel('Time,t', 'FontSize', 12);
ylabel('Height', 'FontSize', 12);

% Drawing Fig.4(c)
figure(2);
len_plot = 1:250;
plot(len_plot, Zt_observ_att(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ2(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ3(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ4(2,len_plot), 'LineWidth', 1.2);
hold on;
a = 160;
b = 210;
x_fill = [a,b,b,a];
y_fill = [min(ylim),min(ylim),max(ylim),max(ylim)];
patch(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
legend('FO via Sensor 1 with attack','FO via Sensor 2','FO via Sensor 3','FO via Sensor 4')
xlabel('Time,t', 'FontSize', 12);
ylabel('Height', 'FontSize', 12);

% Drawing the small figure in Fig.4(c) 
figure(3);
len_plot = 150:220;
plot(len_plot, Zt_observ_att(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ2(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ3(2,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ4(2,len_plot), 'LineWidth', 1.2);
hold on;
a = 160;
b = 210;
x_fill = [a,b,b,a];
y_fill = [min(ylim),min(ylim),max(ylim),max(ylim)];
patch(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time,t', 'FontSize', 12);
ylabel('Height', 'FontSize', 12);

% Drawing Fig.4(d)
figure(4);
len_plot = 1:250;
plot(len_plot, Zt_observ_att(1,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ2(1,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ3(1,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ4(1,len_plot), 'LineWidth', 1.2);
hold on;
a = 160;
b = 210;
x_fill = [a,b,b,a];
y_fill = [min(ylim),min(ylim),max(ylim),max(ylim)];
patch(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
legend('Value of S1 with fault','Recovered value of S1 by S2','Recovered value of S1 by S3','Recovered value of S1 by S4')
xlabel('Time,t', 'FontSize', 12);
ylabel('Height', 'FontSize', 12);

% Drawing the small figure in Fig.4(d) 
figure(5);
len_plot = 150:220;
plot(len_plot, Zt_observ_att(1,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ2(1,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ3(1,len_plot), 'LineWidth', 1.2);
hold on;
plot(len_plot, Zt_observ4(1,len_plot), 'LineWidth', 1.2);
hold on;
a = 160;
b = 210;
x_fill = [a,b,b,a];
y_fill = [min(ylim),min(ylim),max(ylim),max(ylim)];
patch(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time,t', 'FontSize', 12);
ylabel('Height', 'FontSize', 12);

% Parameter Calculation, considering the fluid is oil
function Q = discretized(d,l,dt)                   
rho = 885;                                         % Liquid density (kg/m³)
g = 9.81;                                          % Gravitational acceleration (N/kg)
u = 582.95;                                        % Dynamic viscosity (cp)
Q = dt * rho * g * 3.1415 * d^4 /(128 * l * u);
end


