clear; clc; close all;
addpath('utils')

%% Parameter Settings
sigma = 10;      % Prandtl number
rho = 28;        % Rayleigh number
beta = 8/3;      % Geometric parameter

% Simulation parameters
dt = 0.01;       % Time step
T = 100;         % Total simulation time
num_steps = floor(T/dt);  % Total number of steps

% Initial state
x0 = 1; y0 = 1; z0 = 1;
X = [x0; y0; z0];  % State vector

% Input signal settings
u_amp = 5;         % Input amplitude
u_freq = 0.5;      % Input frequency

% Define output matrix
C = [1, 0, 0];      % Output x
D = zeros(2, 3);     % Direct feedthrough matrix

% Output dimensions
p = size(C, 1);     % Output dimension
m = 3;              % Input dimension
n = 3;              % State dimension

% 3. Preallocate storage
x_history = zeros(n, num_steps+1);
u_history = zeros(m, num_steps);
y_history = zeros(p, num_steps+1);

% Initial values
x_history(:, 1) = X;
y_history(:, 1) = C * X;

%% Main simulation loop
for k = 1:num_steps
    t = (k-1)*dt;
    
    % Generate input signal
    % Can choose different input types
    input_type = 1;  % 1: Sine wave, 2: Step, 3: Random
    
    switch input_type
        case 1
            % Sine wave input
            u = u_amp * [sin(2*pi*u_freq*t);
                        sin(2*pi*u_freq*t + pi/3);
                        sin(2*pi*u_freq*t + 2*pi/3)];
        case 2
            % Step input
            if t < 20
                u = [0; 0; 0];
            else
                u = [1; 0.5; 0.2];
            end
        case 3
            % Random input
            u = 0.1 * randn(3, 1);
        otherwise
            u = zeros(3, 1);
    end
    
    % Store input
    u_history(:, k) = u;
    
    % Current state
    x = x_history(:, k);
    
    % Runge-Kutta 4th order (RK4)
    k1 = dt * lorenz_derivatives(x, sigma, rho, beta, u);
    k2 = dt * lorenz_derivatives(x + k1/2, sigma, rho, beta, u);
    k3 = dt * lorenz_derivatives(x + k2/2, sigma, rho, beta, u);
    k4 = dt * lorenz_derivatives(x + k3, sigma, rho, beta, u);
    
    % Update state
    x_new = x + (k1 + 2*k2 + 2*k3 + k4)/6;
    
    % Store results
    x_history(:, k+1) = x_new;
    y_history(:, k+1) = C * x_new;
end

% Time vector
time = 0:dt:T;

% Save data for observer design
data.x = x_history;
data.y = y_history;
data.u = u_history;
data.time = time;
data.params.sigma = sigma;
data.params.rho = rho;
data.params.beta = beta;
data.params.dt = dt;
data.params.T = T;

fprintf('\nData saved to structure ''data''\n');

% Prepare data for observer
u_d = u_history;         % Input data
y_d = y_history;         % Output data
x_d = x_history(:, 1: end-1);  % State data 

% Target function to estimate 
L = [0, 1, 0];  % L matrix defines z = Lx 
r = size(L, 1);  % Dimension of target function

%% Estimator Design

[Sigma2, flag, newL2] = Reduced_order_observ_augmented_subspace_intersection_v3_partial(u_d, x_d, x_d, C, L, eye(n));

newr2 = size(newL2, 1);

%% Test estimation performance
% Initial state  
x0_test = 0; y0_test = 0; z0_test = 0;
X0_test = [x0_test; y0_test; z0_test];  

% Preallocate test data storage
x_history_test = zeros(3, num_steps+1);
u_history_test = zeros(3, num_steps);
y_history_test = zeros(p, num_steps+1);
z_history_test = zeros(newr2, num_steps);

% Initial values for test set
x_history_test(:, 1) = X0_test;
y_history_test(:, 1) = C * X0_test;
z_history_test(:, 1) = newL2 * X0_test;  % True value of target function

% Generate test data
for k = 1:num_steps
    t = (k-1)*dt;
    
    % Generate input signal
    input_type = 1;
    
    switch input_type
        case 1
            % Sine wave input
            u_test = u_amp * [sin(2*pi*u_freq*t);
                        sin(2*pi*u_freq*t + pi/3);
                        sin(2*pi*u_freq*t + 2*pi/3)];
        case 2
            % Step input
            if t < 20
                u_test = [0; 0; 0];
            else
                u_test = [1; 0.5; 0.2];
            end
        case 3
            % Random input
            u_test = 0.05 * randn(3, 1);
        otherwise
            u_test = zeros(3, 1);
    end
    
    % Store input
    u_history_test(:, k) = u_test;
    
    % Current state
    x_test = x_history_test(:, k);
    
    % Calculate Lorenz derivatives
    % Runge-Kutta 4th order (RK4)
    k1_test = dt * lorenz_derivatives(x_test, sigma, rho, beta, u_test);
    k2_test = dt * lorenz_derivatives(x_test + k1_test/2, sigma, rho, beta, u_test);
    k3_test = dt * lorenz_derivatives(x_test + k2_test/2, sigma, rho, beta, u_test);
    k4_test = dt * lorenz_derivatives(x_test + k3_test, sigma, rho, beta, u_test);
    
    % Update state
    x_new_test = x_test + (k1_test + 2*k2_test + 2*k3_test + k4_test)/6;
    
    % Store results
    x_history_test(:, k+1) = x_new_test;
    y_history_test(:, k+1) = C * x_new_test;
end

% Prepare test data for observer
u_d_test = u_history_test;         % Test input data
x_d_test = x_history_test(:, 1:10000);  % Test state data
y_d_test = y_history_test(:, 1:10000);  % Test output data
ztrue_test = newL2 * x_d_test;

% Observer initialization
zobs_test = zeros(newr2, num_steps); 
zobs_test(:, 1) = rand(newr2, 1);  % Random initial condition

% Error initialization
zerror_test = zeros(newr2, num_steps);
zerror_test(:, 1) = zobs_test(:, 1) - ztrue_test(:, 1);

% Observer simulation
for k = 1:num_steps-1
    Stack(:,k) = [u_d_test(:, k); y_d_test(:, k); y_d_test(:, k+1); zobs_test(:, k)]; 
    zobs_test(:, k+1) = Sigma2 * Stack(:, k);
    zerror_test(:, k+1) = zobs_test(:, k+1) - ztrue_test(:, k+1);
end

% Calculate relative error
for k = 1:num_steps
    for i = 1:r
        z_rerror_test(i, k) = (zobs_test(i, k) - ztrue_test(i, k))/ztrue_test(i, k);
    end
end

%% Visualization
% Figure 1: Comparison of estimated vs true values
plot_steps = 6000;
figure(1)
plot(1:plot_steps, ztrue_test(1, 1:plot_steps), 'b-', 'LineWidth', 2.5)
hold on
plot(1:plot_steps, zobs_test(1, 1:plot_steps), 'r-', 'LineWidth', 2)
xlabel('Time','FontSize',17)
ylabel('State Value','FontSize',17)
legend('True value','Estimated value', 'Location', 'best','FontSize',17)
title('Comparison of True and Estimated States','FontSize',14)
set(gca,'FontSize',16)
grid on

% Figure 2: Phase plane comparison
figure(2)
plot(x_history_test(1,2:end), zobs_test(1,:), 'b-', 'LineWidth', 0.5)
hold on 
plot(x_history_test(1,2:end), x_history_test(2,2:end), 'r-', 'LineWidth', 0.5)
legend('Estimated trajectory (x-y plane)', 'True trajectory (x-y plane)', 'Location', 'best','FontSize',14)
xlabel('x','FontSize',16)
ylabel('y','FontSize',16)
title('Phase Plane Comparison: True vs Estimated','FontSize',14)
set(gca,'FontSize',15)
grid on

%% Lorenz System Derivative Function
function dxdt = lorenz_derivatives(X, sigma, rho, beta, u)
    % Derivative of Lorenz system
    % X: current state [x; y; z]
    % u: input vector [u1; u2; u3]
    
    x = X(1);
    y = X(2);
    z = X(3);
    
    % Lorenz equations
    dxdt = zeros(3,1);
    dxdt(1) = sigma * (y - x) + u(1);      % dx/dt
    dxdt(2) = x * (rho - z) - y + u(2);    % dy/dt
    dxdt(3) = x * y - beta * z + u(3);     % dz/dt
end
