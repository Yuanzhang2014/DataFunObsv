%% DEMO FOR SIMULATION IN FIG.3 - PARTIAL INFORMATION IN SPARSE NETWORKS
% This script provides a demo for numerical simulation in Figure 3 with 
% RRMSE plot example in ER networks with varying partial state information.

%------------- BEGIN CODE --------------

clear 
clc

% Experimental parameters for partial information study
m = 80;                          % Number of inputs
p = 100;                         % Number of outputs
r = 10;                          % Number of functionals to estimate
n = 250;                         % System dimension
num_nodes = 9;                   % Number of partial information levels
step = round(2 * n / (num_nodes - 1));  % Step size for information levels
num_tries = 100;                 % Monte Carlo simulations
len_data = 800;                  % Training data length
p_diag = 0.6;                    % Diagonal probability for system generation
p_threshold = [0.2 0.4 0.6 0.8]; % Edge probability thresholds for ER networks

% Information levels: number of measured states (from 0 to 2n)
W_lib = [0:step:2*n 2*n];

% Initialize result arrays
flags = zeros(length(W_lib), num_tries, length(p_threshold));      % Design success flags
partial_ord = zeros(length(W_lib), num_tries, length(p_threshold)); % Observer orders
RRMSE = zeros(length(W_lib), num_tries, length(p_threshold));      % Performance metrics
time = zeros(length(W_lib), num_tries, length(p_threshold));       % Computation times

% Results file with timestamp
str_save = ['partial_info_sparse_ER_' num2str(num_tries) 'tries_' datestr(now, 'mm_dd_yy') '.mat'];
str_save = ['./' str_save];

%% Main simulation loop
t_start = datetime('now');

for w_idx = 1:length(W_lib)
    Wlen = W_lib(w_idx);  % Number of measured states in this iteration
    
    for p_idx = 1:length(p_threshold)
        for count = 1:num_tries
            % Generate observable ER system with given edge probability
            while 1
                [A, B, C, L, ~, IsObsv] = sparse_ER_gen(n, m, p, r, p_threshold(p_idx), p_diag);
                if IsObsv
                    break;  % Ensure we only use observable systems
                end
            end
            
            % Generate training trajectory
            [U_train, X_train] = model2traj(A, B, len_data, 'train');
            
            % Create partial state measurement matrix
            W_idx = randperm(2 * n, Wlen);      % Random selection of measured states
            identity = eye(2 * n);
            W_transform = identity(W_idx, :);   % Measurement transformation
            W = X_train(W_idx, :);              % Partial state measurements
            
            % Reduced-order observer design with partial information
            tic
            [Sigma2, flag, newL2] = Reduced_order_observ_augmented_subspace_intersection_v3_partial(U_train, X_train, W, C, L, W_transform);
            time(w_idx, count, p_idx) = toc;    % Record computation time
            flags(w_idx, count, p_idx) = flag;  % Store design success flag
            partial_ord(w_idx, count, p_idx) = size(Sigma2, 1);  % Store observer order
            
            % Performance evaluation on test data
            [U_test, X_test] = model2traj(A, B, 2 * len_data, 'test');
            Z_test = L * X_test;  % True functional values
            
            % Compute RRMSE (Relative Root Mean Square Error)
            error_data = Simulation_agumentation_observer_with_trajectory_no_plot(Sigma2, X_test, U_test, C, newL2, L);
            acc_error = sqrt(trace(error_data(1:r, end-800+1:end)' * error_data(1:r, end-800+1:end)) / ...
                        trace(Z_test(1:r, end-800+1:end)' * Z_test(1:r, end-800+1:end)));
            RRMSE(w_idx, count, p_idx) = acc_error;
            
            % Display current result
            disp('accumulative error:');
            disp(acc_error);
            
            % Periodic saving to preserve progress
            if mod(count, 10) == 0
                save(str_save);
            end
        end
    end
end

t_end = datetime('now');
total_time = between(t_start, t_end);

disp('total time spent:')
disp(total_time)

%% Performance visualization
% Process results and generate RRMSE plot

% Handle numerical instability cases (Inf/NaN values)
RRMSE(isnan(RRMSE)) = 0;
mask = RRMSE > 0;
e_avg = sum(RRMSE .* mask, 2) ./ sum(mask, 2);
e_avg = reshape(e_avg, size(e_avg, 1), size(e_avg, 3));

% Color scheme for different edge probability thresholds
colors = [255, 204, 51;
          27, 158, 119;
          166, 86, 40;
          117, 112, 179;
          231, 41, 138] / 255;

%%
load('ER_error.mat')
% Create RRMSE vs partial information plot
figure
x_plot = [5:round((100 - 5) / 8):100 100];  % X-axis as percentage of measured states

for i = 1:4
    p = plot(x_plot, e_avg(:, i), 'LineWidth', 1, 'Color', colors(5 - i, :));
    hold on
end

% Logarithmic scale for better visualization of error range
set(gca, 'Yscale', 'log')

% Figure formatting
set(gca, 'Box', 'on', ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'in', 'TickLength', [.012 .012], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [0 0 0], 'YColor', [0 0 0], ...
         'FontName', 'Arial', ...
         'FontSize', 10)

% Axis limits and ticks
xl = xlim;
xlim([xl(1) - 2, xl(2) + 2]);
set(gca, 'xtick', 5:45:95);
set(gca, 'ytick', [10^(-5), 10^0, 10^5]);
ylim([6 * 10^(-6), 5 * 10^6]);
set(gcf, 'Position', [200 200 360 280]);

% Final plot styling
ax = gca;
ax.YColor = [0, 0, 0];
ax.XColor = [0, 0, 0];
ax.LineWidth = 0.6;


%------------- END CODE --------------
