%% Target estimation in nonlinear power grids
% 39-node system (10 generator nodes and 29 load nodes)
% Other m-files required: Reduced_order_observ_augmented_subspace_intersection_affine_v3
% Darouch_observer_with_trajectory_affine_v1.m
% IsObservable.m
% rand_svd.m
% data_functional_observability_test.m
% power39PB.m  
% Data presents network information for IEEE 39, with topology data sourced from the following reference and values randomly generated.
% Gérin-Lajoie, L., Saad, O., & Mahseredjian, J. (2015). IEEE PES Task Force on Benchmark Systems for Stability Controls [R]. EMTP-RV 39-bus system.


                     
% Author: Wenxuan Xu
% email:  xuwenxuan1209@163.com
% Last revision: Oct-28-2025


clear; close all; clc;
addpath('../../utils')
%% Parameter settings
rng(30); % Set random seed for reproducibility

% System parameters
n_g = 10;      % Number of generator nodes
n_l = 29;      % Number of load nodes
N = n_g + n_l; % Total number of nodes
dt = 0.001;    % Time step
T = 5;         % Total simulation time
steps = T/dt;  % Total number of simulation steps

% Network topology parameters
load power39PB.mat

J = [zeros(n_g,n_l), eye(n_g); eye(n_l), zeros(n_l,n_g)];
G = J * B * inv(J);

V = ones(N,1); % Voltage magnitudes 
K = G;
for i = 1:N
    for j = 1:N
        K(i,j) = V(i)*V(j)*G(i,j);
    end
end

% Generator parameters (second-order oscillators)
M = zeros(N,1); % Inertia constants
D = zeros(N,1); % Damping constant
P = zeros(N,1); % Power injection

% Generator node parameters (nodes 1–10), randomly generated within a reasonable range
M(1:n_g) = 1 + 0.5*rand(n_g,1);   % Inertia constants in [1, 1.5]
D(1:n_g) = 0.5 + 0.3*rand(n_g,1); % Damping constant in [0.5, 0.8]
P(1:n_g) = 0.8 + 0.4*rand(n_g,1); % Power injection in [0.8, 1.2]

% Load node parameters (nodes 11–39)
D(n_g+1:end) = 0.3 + 0.2*rand(n_l,1);   % Damping constant in [0.3, 0.5]
P(n_g+1:end) = -0.5 - 0.3*rand(n_l,1);  % Power injection in [-0.8, -0.5]

% Power generation
P_total = sum(P);
if P_total > 0
    P(1:n_g) = P(1:n_g) - P_total/n_g;
else
    P(1:n_g) = P(1:n_g) - P_total/n_g;
end

%% Initialization of state variables
phi = zeros(N, steps);       % Phase angles
omega = zeros(N, steps);     % Angular frequencies (only for generator nodes)
dphi_dt = zeros(N, steps);   % Phase derivatives
domega_dt = zeros(N, steps); % Angular frequency derivatives

% Random initial conditions
pi0 = 2*pi*rand(N,1) - pi;
phi(:,1) = pi0;              % Initial phases in [-π, π]
omega0 = 0.1*randn(n_g,1);
omega(1:n_g,1) = omega0;     % Initial generator angular velocities

%% Generate simulation data
u_traj = [];
for t = 1:steps-1
    % Random input
    u = 0.01*rand(N,1);
    
    % Compute coupling term at the current time step
    coupling = zeros(N,1);
    for i = 1:N
        for j = 1:N
            if i ~= j && K(i,j) ~= 0
                coupling(i) = coupling(i) + K(i,j) * sin(phi(j,t) - phi(i,t));
            end
        end
    end
    
    % Update generator nodes
    for i = 1:n_g
        domega_dt(i,t) = (P(i) - D(i)*omega(i,t) + coupling(i)) / M(i);
        dphi_dt(i,t) = omega(i,t);
        
        % Forward Euler discretization
        omega(i,t+1) = omega(i,t) + dt * domega_dt(i,t);
        phi(i,t+1) = phi(i,t) + dt * (dphi_dt(i,t) + u(i));
    end
    
    % Update load nodes
    for i = n_g+1:N
        dphi_dt(i,t) = (P(i) + coupling(i)) / D(i);
        phi(i,t+1) = phi(i,t) + dt * (dphi_dt(i,t) + u(i));
        omega(i,t+1) = 0;  % Load nodes have no frequency dynamics
    end
    
    u_traj = [u_traj, u];
end

x_traj = [omega(1:n_g,:); phi(1:end,:)];
x_traj = x_traj(:,1:end-1);

%% Outputs and measurement nodes
p = 20; r = 2;
I = eye(49);
crows = [11,14,16,17,18,19,20,21,25,27,28,30,31,35,37,38,40,41,44,49];
C = I(crows,:);
lrows = [2,4,13,15,32,34];
L = I(lrows,:);

%% Data-driven functional observer design
[Sigma, exist, newL] = Reduced_order_observ_augmented_subspace_intersection_affine_v3(u_traj, x_traj, C, L);

%% Simulation validation
phi = zeros(N, steps);     
omega = zeros(N, steps);   
dphi_dt = zeros(N, steps);  
domega_dt = zeros(N, steps); 

% Initial state for testing
phi(:,1) = pi0 + 0.01*randn(N,1);  
omega(1:n_g,1) = omega0 + 0.01*randn(n_g,1); 

u_test = [];
for t = 1:steps-1
    u = 0.01*rand(N,1);
    coupling = zeros(N,1);

    for i = 1:N
        for j = 1:N
            if i ~= j && K(i,j) ~= 0
                coupling(i) = coupling(i) + K(i,j) * sin(phi(j,t) - phi(i,t));
            end
        end
    end
    
    for i = 1:n_g
        domega_dt(i,t) = (P(i) - D(i)*omega(i,t) + coupling(i)) / M(i);
        dphi_dt(i,t) = omega(i,t);
        omega(i,t+1) = omega(i,t) + dt * domega_dt(i,t);
        phi(i,t+1) = phi(i,t) + dt*(dphi_dt(i,t) + u(i));
    end
    
    for i = n_g+1:N
        dphi_dt(i,t) = (P(i) + coupling(i)) / D(i);
        phi(i,t+1) = phi(i,t) + dt*(dphi_dt(i,t) + u(i));
        omega(i,t+1) = 0;
    end
    
    u_test = [u_test, u];
end

x_test = [omega(1:n_g,:); phi(1:end,:)];
x_test = x_test(:,1:end-1);

% simulation
X = x_test;
X0 = X(:,1);
Y0 = C*X0;
r = size(L,1);
observer_Z0 = 10*rand(r,1);
observer_Zall = observer_Z0;
Zall = L*X;

for t = 2:size(u_traj,2)
    U = u_traj(:,t-1);
    Y1 = C*X(:,t);
    observer_Z1 = Sigma * [U; Y0; Y1; observer_Z0; ones(1,size(U,2))];
    observer_Z0 = observer_Z1;
    Y0 = Y1;
    observer_Zall = [observer_Zall observer_Z1];
end

Z_error = Zall - observer_Zall;

%% plot results
t = (1:size(Zall,2)) * 0.001;

% First figure: comparison of row 1 data
figure('Name','Comparison of Row 1 Data');
plot(t, Zall(1,:), 'b-', 'LineWidth',1.2); hold on;
plot(t, observer_Zall(1,:), 'r--', 'LineWidth',1.2); hold on;
plot(t, Zall(2,:), 'g-', 'LineWidth',1.2); hold on;
plot(t, observer_Zall(2,:), 'y--', 'LineWidth',1.2);
legend({'$\dot{\phi}_2(t)$', '$\dot{\hat{\phi}}_2(t)$', '$\dot{\phi}_4(t)$', '$\dot{\hat{\phi}}_4(t)$'}, 'Interpreter','latex');

% Second figure: comparison of row 2 data
figure('Name','Comparison of Row 2 Data');
plot(t, Zall(3,:), 'b-', 'LineWidth',1.2); hold on;
plot(t, observer_Zall(3,:), 'r--', 'LineWidth',1.2); hold on;
plot(t, Zall(4,:), 'g-', 'LineWidth',1.2); hold on;
plot(t, observer_Zall(4,:), 'y--', 'LineWidth',1.2);
legend({'${\phi}_3(t)$', '$\hat{\phi}_3(t)$', '${\phi}_5(t)$', '$\hat{\phi}_5(t)$'}, 'Interpreter','latex');

% Third figure: comparison of row 3 data
figure('Name','Comparison of Row 3 Data');
plot(t, Zall(5,:), 'b-', 'LineWidth',1.2); hold on;
plot(t, observer_Zall(5,:), 'r--', 'LineWidth',1.2); hold on;
plot(t, Zall(6,:), 'g-', 'LineWidth',1.2); hold on;
plot(t, observer_Zall(6,:), 'y--', 'LineWidth',1.2);
legend({'${\phi}_{22}(t)$', '$\hat{\phi}_{22}(t)$', '${\phi}_{24}(t)$', '$\hat{\phi}_{24}(t)$'}, 'Interpreter','latex');
% plot 4: Network topology
figure
G = graph(K);
p = plot(G, 'Layout','force', 'UseGravity',true);
p.NodeCData = [ones(n_g,1); 2*ones(n_l,1)]; % Generator nodes in red, load nodes in blue
p.MarkerSize = 4;
colormap([1 0 0; 0 0 1]); % Red: generators, Blue: loads
