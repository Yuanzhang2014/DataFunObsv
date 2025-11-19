function [Z_error, R_t_sigma_o] = Simulation_agumentation_observer_with_trajectory_pure_data(Sigma, x_traj, y_traj, u_traj, z_traj)
% =========================================================================
% Filename: Simulation_agumentation_observer_with_trajectory_pure_data.m
% Purpose: Simulate and evaluate the performance of a data-driven augmented state observer
% Inputs: 
%   Sigma - Observer parameter matrix
%   x_traj - State trajectory
%   y_traj - Output trajectory  
%   u_traj - Input trajectory
%   z_traj - Target functional state trajectory
% Outputs:
%   Z_error - Observer error trajectory
%   R_t_sigma_o - Root mean square error performance metric
% Main functions:
%   1. Simulate observer dynamics using Darouch observer parameters
%   2. Calculate observer errors and performance metrics
%   3. Generate multiple visualization plots for performance analysis
% =========================================================================

X = x_traj;
Y = y_traj;
X0 = X(:, 1);
Y0 = Y(:, 1); 
r0 = size(z_traj, 1);          
r = size(Sigma, 1);            


observer_Z0 = rand(r, 1);       
observer_Z = observer_Z0(1:r0, :); 
Zall = z_traj;               

% Observer simulation loop
for t = 2:size(u_traj, 2)
    U = u_traj(:, t-1);         
    Y1 = y_traj(:, t);          
    observer_Z1 = Sigma * [U; Y0; Y1; observer_Z0];
    observer_Z0 = observer_Z1;
    Y0 = Y1;
    observer_Z = [observer_Z, observer_Z1(1:r0, :)];
end

Z_error = Zall - observer_Z; 
R_t_sigma_o = sqrt(sum(Z_error.^2, 1) / (r0^2)); 

%% Visualization section - Generate multiple analysis plots
Len_data = size(z_traj, 2) - 1;
Len_data_cut = 1;
e_num_rows = r0;

% Figure 1: Absolute error trajectories
figure(1)
for i = 1:e_num_rows
    subplot(e_num_rows, 1, i); 
    plot(Len_data_cut:Len_data, Z_error(i, Len_data_cut:Len_data), 'LineWidth', 1);
    xlabel('{Time}', 'Fontsize', 12);
    ylabel({['$e_{', num2str(i), '}(t)$']}, 'Interpreter', 'latex', 'Fontsize', 12);

    if i == 1
        title('Errors in a data-driven functional observer', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'b');
    end
    hold on;
end

% Figure 2: Relative error trajectories  
figure(2)
rZ_error = Z_error ./ Zall;
for i = 1:e_num_rows
    subplot(e_num_rows, 1, i); 
    plot(1:Len_data, rZ_error(i, 1:Len_data), 'LineWidth', 1);
    xlabel('{Time}', 'Fontsize', 16);
    ylabel({['$\delta e_{', num2str(i), '}(t)$']}, 'Interpreter', 'latex', 'Fontsize', 16);
    
    if i == 1
        title('Relative errors in a data-driven functional observer', 'FontWeight', 'bold', 'FontSize', 10, 'Color', 'b');
    end
    hold on;
end

% Figure 3: Comparison of estimated vs true states (subplot format)
figure(3)
for i = 1:e_num_rows
    subplot(e_num_rows, 1, i); 
    plot(1:Len_data, observer_Z(i, 1:Len_data), 'o', 'LineWidth', 1); 
    hold on;
    plot(1:Len_data, Zall(i, 1:Len_data), '--', 'LineWidth', 1);   
    hold on;
    xlabel('Time Step', 'Fontsize', 16);
    legend({['$\hat z_{', num2str(i), '}(t)$'], ['$z_{', num2str(i), '}(t)$']}, 'Interpreter', 'latex', 'Fontsize', 16);
end

% Figure 4: Overlaid trajectory comparison with color coding
figure(4)
c = colormap(jet(e_num_rows)); 
for k = 1:e_num_rows
    plot(1:Len_data, Zall(k, 1:Len_data), '--', 'color', c(k, :), 'LineWidth', 1);  % True states
    hold on;
    plot(1:Len_data, observer_Z(k, 1:Len_data), '-', 'color', c(k, :), 'LineWidth', 2);  % Estimated states
    hold on;
end
xlabel('Time ', 'Fontsize', 16);
leg_str = cell(1, 2*e_num_rows);
for k = 1:e_num_rows
    leg_str{2*k-1} = ['$z_{', num2str(k), '}(t)$'];   
    leg_str{2*k} = ['$\hat z_{', num2str(k), '}(t)$']; 
end
legend(leg_str, 'Interpreter', 'latex', 'Location', 'best', 'Fontsize', 16);