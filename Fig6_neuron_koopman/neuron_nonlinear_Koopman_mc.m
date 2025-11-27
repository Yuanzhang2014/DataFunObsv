% =========================================================================
% Filename: neuron_nonlinear_Koopman_mc.m
% Purpose: Main script for comparing observer performance across different Koopman lifting dimensions
% Inputs: None (all parameters defined internally)
% Outputs: Performance comparison plots for different lifting dimensions
% Main functions:
%   1. Generate training and test data for neuron network system
%   2. Apply Koopman operator lifting with different dimensions
%   3. Design and test data-driven observers
%   4. Compare observer performance across different lifting dimensions
% =========================================================================

clear;
clc;
close all;
rng(0);
addpath('../utils')

%% Generate training dataset
node_num = 50;                 
p = 80;                        
r = 5;                      
p_edge = 0.3;                 

% Define different lifting dimensions
lift_dims = [30, 50, 100, 200]; 
num_monte_carlo = 100; % Number of Monte Carlo experiments

% Generate ER random network topology
G = rand(node_num) < p_edge;
Topology = double(G);

n = node_num * 3;             
m = node_num;                   
N = n * m + 2 * n + 1;       
I = eye(n);

% Randomly select rows for C and L matrices
rows = randperm(n, p+r);
crows = rows(:, 1:p);   
lrows = rows(:, p+1:p+r);   
C = I(crows, :);               
L = I(lrows, :);               

X0 = 10 * rand(n, 1);          
% Generate nonlinear system data
[X_non, Y_non, Z_non, U] = neuron_nonlinear_centerpoint(node_num, p, r, Topology, C, L, N, X0);

%% Koopman lifting process
for idx = 1:length(lift_dims)
    num_centers_x = lift_dims(idx);  % Current lifting dimension
    num_centers_y = 0;               % No lifting in Y direction
    centers_x = 1.5 * lhsdesign(num_centers_x, 2) - 0.75;
    centers_y = 1.5 * lhsdesign(num_centers_y, 2) - 0.75;
    num_steps = N;
    
    % Perform Koopman lifting on training data
    [Phi_x, Phi_y] = lift_koopman(X_non, Y_non, num_centers_x, num_centers_y, num_steps, centers_x, centers_y);
    
    %% Observer design section
    % Prepare trajectory data
    x_traj = Phi_x; 
    y_traj = Phi_y; 
    z_traj = Z_non;
    u_traj = U;
    Up = u_traj;
    
    % Design Darouch observer using pure data
    [Sigma, order] = Darouch_observer_with_trajectory_pure_data(u_traj, y_traj, z_traj);
    
    % Prepare data matrices for observer verification
    up = Up(:, 1:end-1);
    xp = x_traj(:, 1:end-1);
    xf = x_traj(:, 2:end);
    yp = y_traj(:, 1:end-1);
    yf = y_traj(:, 2:end);
    zp = xp;
    zf = xf;
        
    % Observer design matrix construction
    P = [up; yp; yf; zp];    
    f = [P; zf];
    M1 = P;

    % Check persistent excitation condition
    if rank(M1) < size(M1, 1)
        warning('Data does not satisfy persistent excitation condition');
    end
    
    % Check system consistency
    check_consistent(P, xf)     % Verify rank_svd(P) == rank_svd(f)
    
    %% Monte Carlo experiments with different initial conditions
    fprintf('Starting Monte Carlo experiments for lifting dimension %d...\n', lift_dims(idx));
    
    % Preallocate storage for error trajectories
    error_trajectories = zeros(num_monte_carlo, N);
    successful_runs = 0;
    
    for N_r = 1:num_monte_carlo
        try
            % Generate different initial condition for each experiment
            probability = 0.9;     
            noiseMask = (rand(size(X0)) < probability); 
            noise_intensity = 0.01 + 0.02 * rand(); % Vary noise intensity slightly
            noise = noise_intensity * randn(size(X0)) .* noiseMask; 
            noisyX0 = X0 + noise;  

            [X_non_test, Y_non_test, Z_non_test, U_test] = neuron_nonlinear_centerpoint(node_num, p, r, Topology, C, L, N, noisyX0);
            
            [Phi_x_test, Phi_y_test] = lift_koopman(X_non_test, Y_non_test, num_centers_x, num_centers_y, num_steps, centers_x, centers_y);

            x_traj_test = Phi_x_test;
            y_traj_test = Phi_y_test;
            u_traj_test = U_test;
            z_traj_test = Z_non_test;
            
            [z_error, R_t_sigma_o] = Simulation_agumentation_observer_with_trajectory_pure_data(Sigma, x_traj_test, y_traj_test, u_traj_test, z_traj_test);
            
            % Store the error trajectory for this experiment
            if length(R_t_sigma_o) >= N
                error_trajectories(N_r, :) = R_t_sigma_o(1:N);
            else
                % Pad with last value if trajectory is shorter
                error_trajectories(N_r, :) = [R_t_sigma_o, R_t_sigma_o(end)*ones(1, N-length(R_t_sigma_o))];
            end
            
            successful_runs = successful_runs + 1;
            
        catch ME
            fprintf('Experiment %d failed: %s\n', N_r, ME.message);
            error_trajectories(N_r, :) = NaN;
        end
    end
    
    fprintf('Completed %d successful experiments for lifting dimension %d.\n', successful_runs, lift_dims(idx));
    
    %% Calculate RMSE according to the formula: R_t = sqrt(1/N_r * sum(||x_t^i - x_hat_t^i||^2))
    % Remove failed experiments
    valid_experiments = error_trajectories(~all(isnan(error_trajectories), 2), :);
    num_valid_runs = size(valid_experiments, 1);
    
    if num_valid_runs == 0
        error('No successful experiments for lifting dimension %d.', lift_dims(idx));
    end
    
    % Calculate RMSE according to the formula
    % R_t = sqrt( (1/N_r) * sum_{i=1}^{N_r} ||x_t^i - x_hat_t^i||^2 )
    squared_errors = valid_experiments.^2;  % Square the errors
    mean_squared_errors = mean(squared_errors, 1);  % Average over experiments
    RMSE_trajectory = sqrt(mean_squared_errors);  % Square root to get RMSE
    
    % Store the RMSE trajectory for this lifting dimension
    R_t_cell{idx} = RMSE_trajectory;
    
    % Also store individual error trajectories for envelope calculation
    error_cell{idx} = valid_experiments;
end

%% Performance visualization - RMSE with envelope
N_test = 6000;

% Figure 1: RMSE comparison
figure('Position', [100, 100, 800, 600]);
hold on;
colors = lines(length(lift_dims));    
line_styles = {'-', '--', ':'};
for idx = 1:length(lift_dims)
    RMSE_traj = R_t_cell{idx};
    time_steps = 1:N_test;
    plot(time_steps, RMSE_traj(1:N_test), ...
         'Color', colors(idx, :), ...
         'LineStyle', line_styles{mod(idx-1, length(line_styles)) + 1}, ...
         'LineWidth', 2, ...
         'DisplayName', ['Lift Dim = ', num2str(lift_dims(idx))]);
end
hold off;
xlabel('Time Step', 'FontSize', 12);
ylabel('RMSE', 'FontSize', 12);
%title('RMSE Comparison (Averaged over Experiments)', 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'best');
grid on;


