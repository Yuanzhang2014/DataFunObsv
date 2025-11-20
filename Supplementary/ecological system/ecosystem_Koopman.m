% Target estimation in ecological systems.
% Example: An ecological system with 20 interacting species.
% Other m-files required: Darouch_observer_with_trajectory_pure_data.m
%                         lift_koopman.m
%                         rank_svd.m
%                         Simulation_agumentation_observer_with_trajectory_pure_data_ecos.m
% Random generation of data parameters

                     
% Author: Wenxuan Xu
% email:  xuwenxuan1209@163.com
% Last revision: Oct-28-2025

clear;
clc;
close all;
rng(0);
addpath('../../utils')
% Define different lifting dimensions 
lift_dims = [100];
R_t_cell = cell(length(lift_dims), 1);  % Store R_t_sigma_o results
n = 20;   % Number of species
m = 20;   % Dimension of input

% Intrinsic growth rate vector R (per-species intrinsic growth/decay)
R = [1.5;    % Stipa grandis
     1.2;    % Leymus chinensis
     2.0;    % Setaria viridis
     0.9;    % Medicago sativa
     0.3;    % Caragana
     0.4;    % Stellera chamaejasme
     0.7;    % Artemisia spp
     1.0;    % Cyperaceae
     0.8;    % Grasshopper
     0.6;    % Meadow vole
     0.4;    % Plateau pika
     0.2;    % Mongolian gazelle
     0.3;    % Sheep
     0.5;    % Zokor
    -0.1;    % Buzzard
    -0.2;    % Sand fox
    -0.3;    % Wolf
    -0.4;    % Tick
     0.1;    % Dung beetle
     1.2];   % Nematode

% Interaction matrix A (Lotka–Volterra inter-species interactions)
A = zeros(20, 20);

% Set diagonal elements (intra-species competition, must be negative)
diagVals = [-0.025; -0.020; -0.050; -0.015; -0.006; -0.010; -0.014; -0.012;
            -0.100; -0.150; -0.080; -0.040; -0.060; -0.200; -0.008; -0.010;
            -0.005; -0.020; -0.050; -0.300];
A(1:20+1:end) = diagVals;

% Set off-diagonal interaction terms
A(1,2) = -0.005;  % L. chinensis suppresses S. grandis
A(2,1) = -0.003;  % S. grandis suppresses L. chinensis

% Grasshopper (9) feeds on plants (1–3)
A(1,9) = -0.010;  A(9,1) = +0.010;
A(2,9) = -0.010;  A(9,2) = +0.010;
A(3,9) = -0.030;  A(9,3) = +0.010;

% Herbivores eating plants
A(1,11) = -0.025; A(11,1) = +0.008;
A(2,11) = -0.025; A(11,2) = +0.008;

% Predator–prey between mammals
A(10,16) = -0.100; A(16,10) = +0.015; % Vole–fox
A(12,17) = -0.150; A(17,12) = +0.012; % Gazelle–wolf

% Parasitism
A(12,18) = -0.080; A(18,12) = +0.020; % Gazelle–tick

% Decomposer–plant mutualism
A(1,19) = +0.002;  % Dung beetle improves soil nutrient for S. grandis
A(2,19) = +0.003;
A(1,20) = +0.002;  % Nematode improves soil for S. grandis
A(2,20) = +0.003;
A(3,20) = +0.001;

A(10,11) = -0.020; % Pika suppresses vole
A(11,10) = -0.015; % Vole suppresses pika

A(10,15) = -0.100; A(15,10) = +0.018; % Vole–raptor

%% Generate training dataset using nonlinear Lotka–Volterra dynamics
dt = 0.01;
t_end = 50;
T = t_end / dt;
N = max(n * m + 2 * n + 1, T); % Ensure long enough trajectory length

x = zeros(n, N);
u_traj = [];

% Random initial condition
x0 = 200 * rand(n,1);
x(:,1) = x0;

% Forward Euler discretization of dynamics
for k = 1:N
    x_current = x(:,k);
    u = 0.1 * rand(n,1);                            % Random input
    growth_term = x_current .* (R + A * x_current); 
    x(:,k+1) = x_current + dt * (growth_term + u);
    u_traj = [u_traj, u];
end

% training data
X_non = x(:,1:N);
U = u_traj;

% Define output and target matrices (measured species & target species)
r = 2;
I = eye(n);
C = [I(1:10,:); I(12:16,:); I(18:n,:)]; % Measured species
L = [I(11,:); I(17,:)];                 % Target species to estimate

Y_non = C * X_non;
Z_non = L * X_non;

% Koopman lifting
for idx = 1:length(lift_dims)
    num_centers_x = lift_dims(idx);  % Lifting dimension for X
    num_centers_y = 0;               % No lifting for Y

    num_steps = N;
    [Phi_x, Phi_y, centers_x, centers_y] = lift_koopman( ...
        X_non, Y_non, num_centers_x, num_centers_y, num_steps);
    x_traj = Phi_x;
    y_traj = Phi_y;
    z_traj = Z_non;
    u_traj = U;
  % Data-driven Darouach functional observer 
    [Sigma, order] = Darouch_observer_with_trajectory_pure_data( ...
        u_traj, y_traj, z_traj);

  % Generate test data
    x = zeros(n, N);
    x(:,1) = x0 + 40 * (-0.5*ones(n,1) + rand(n,1));  % Perturbed initial state
    u_traj = [];

    for k = 1:N
        x_current = x(:,k);
        u = 0.1 * rand(n,1);
        growth_term = x_current .* (R + A * x_current);
        x(:,k+1) = x_current + dt * (growth_term + u);
        u_traj = [u_traj, u];
    end

    X_non_test = x(:,1:end-1);
    Y_non_test = C * X_non_test;
    Z_non_test = L * X_non_test;
    U_test = u_traj;

    % Apply same lifting to test data
    [Phi_x_test, Phi_y_test] = lift_koopman( ...
        X_non_test, Y_non_test, num_centers_x, num_centers_y, num_steps, centers_x, centers_y);

    % Prepare test trajectories
    x_traj_test = Phi_x_test;
    y_traj_test = Phi_y_test;
    u_traj_test = U_test;
    z_traj_test = Z_non_test;

    % Simulate observer on test data
    [z_error, R_t_sigma_o] = Simulation_agumentation_observer_with_trajectory_pure_data_ecos( ...
    Sigma, x_traj_test, y_traj_test, u_traj_test, z_traj_test);

    R_t_cell{idx} = R_t_sigma_o;
end
%%
% Compute RMSE over 100 independent experiments
% Statistical error evaluation
Rs = [];
for i = 1:100
    x = zeros(n, N);

    xtrain0 = x0;
    % Random initial condition
    x(:,1) = xtrain0 + 10 * (-0.5*ones(n,1) + 1.*rand(n,1));  
    
    u_traj = [];
    for k = 1:N
        x_current = x(:,k);
        u = 0.1 * rand(n,1);
        growth_term = x_current .* (R + A * x_current);
        x(:,k+1) = x_current + dt * (growth_term + u);
        u_traj = [u_traj, u];
    end

    X_non_test = x(:,1:end-1);
    Y_non_test = C * X_non_test;
    Z_non_test = L * X_non_test;
    U_test = u_traj;

    % Koopman lifting for test data
    [Phi_x_test, Phi_y_test] = lift_koopman( ...
        X_non_test, Y_non_test, num_centers_x, num_centers_y, ...
        N, centers_x, centers_y);

    x_traj_test = Phi_x_test;
    y_traj_test = Phi_y_test;
    u_traj_test = U_test;
    z_traj_test = Z_non_test;

    % Observer simulation for test trajectories
    [z_error, R_t_sigma_o] = ...
        Simulation_agumentation_observer_with_trajectory_pure_data_ecos( ...
        Sigma, x_traj_test, y_traj_test, u_traj_test, z_traj_test);

    Rs = [Rs; R_t_sigma_o];

    % Remove rows containing NaNs for valid statistics
    RRs = Rs(~any(isnan(Rs), 2), :);
end

%% Plot 
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultTextFontSize', 10);
set(0, 'DefaultLegendFontSize', 10);
figure(100);
hold on;

meanVals = mean(RRs);

% Compute standard error (SE)
n1 = size(RRs, 1);
stderr = std(RRs) / sqrt(n1);

% 95% confidence interval
ci95 = 1.96 * stderr;

% x-axis indices
x = 1:size(RRs,2);

% Plot 95% CI as filled area
fill([x fliplr(x)], [meanVals+ci95 fliplr(meanVals-ci95)], ...
    [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Plot mean curve
plot(x, meanVals, 'r-', 'LineWidth', 0.5);

set(gcf, 'Position', [200 200 250 200]); 
box on; grid on;
