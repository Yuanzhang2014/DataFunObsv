function [Z_error, R_t_sigma_o] = Simulation_agumentation_observer_with_trajectory_pure_data_ecos(Sigma, x_traj, y_traj, u_traj, z_traj)
%%% Generate observer error using the Darouch observer parameters
%%% It is assumed that L corresponds to the first several rows of Sigma
% Author: Wenxuan Xu
% email:  xuwenxuan1209@163.com
% Last revision: Oct-28-2025

%------------- BEGIN CODE --------------
X  = x_traj;
Y  = y_traj;
X0 = X(:, 1);
Y0 = Y(:, 1);

r0 = size(z_traj, 1);   % Dimension of the target states
r  = size(Sigma, 1);    % Dimension of the augmented states

observer_Z0 = 0 * rand(r,1);     % Initial augmented observer state
observer_Z  = observer_Z0(1:r0); % Initial estimated target state
Zall = z_traj;                   % Ground-truth target trajectories

for t = 2:size(u_traj, 2)
    U  = u_traj(:, t-1);
    Y1 = y_traj(:, t);

    % Update augmented observer state
    observer_Z1 = Sigma * [U; Y0; Y1; observer_Z0];
    observer_Z0 = observer_Z1;
    Y0 = Y1;

    % estimated target state
    observer_Z = [observer_Z observer_Z1(1:r0, :)];
end

% Estimation error
Z_error = Zall - observer_Z;

% (RMSE)
R_t_sigma_o = sqrt(sum(Z_error.^2, 1) / r0);

%%% Plotting error 

Len_data = size(z_traj, 2) - 1;
Len_data = 2000;

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultTextFontSize', 10);
set(0, 'DefaultLegendFontSize', 10);

figure(1)
e_num_rows = r0;

% Plot error time series for each target state
for i = 1:e_num_rows
    plot(1:Len_data, Z_error(i, 1:Len_data), 'LineWidth', 1);
    set(gcf, 'Position', [200 200 250 200]);
    legend({['$e_{1}(t)$'], ['$e_{2}(t)$']}, ...
        'Interpreter', 'latex', 'FontSize', 14);
    hold on;
end

figure(4)
c = colormap(jet(e_num_rows)); 

% Plot true and estimated trajectories
for k = 1:e_num_rows
    plot(1:Len_data, Zall(k, 1:Len_data), '--', ...
         'color', c(k,:), 'LineWidth', 1);   % True trajectories
    hold on

    plot(1:Len_data, observer_Z(k, 1:Len_data), '-', ...
         'color', c(k,:), 'LineWidth', 1);   % Estimated trajectories
    hold on;
end

% Legend construction
leg_str = cell(1, 2 * e_num_rows);  
for k = 1:e_num_rows
    leg_str{2*k - 1} = ['$z_{', num2str(k), '}(t)$'];          % True state
    leg_str{2*k}     = ['$\hat z_{', num2str(k), '}(t)$'];     % Estimated state
end

legend(leg_str, 'Interpreter', 'latex', 'Location', 'best', 'FontSize', 14); 
set(gcf, 'Position', [200 200 250 200]);

end
