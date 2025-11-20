function Z_error = pure_data_simulation_no_plot(sigma, u_traj, y_traj, z_traj)
% PURE_DATA_SIMULATION_NO_PLOT - Simulate functional observer performance
% Simulates the functional observer using the coefficient matrix sigma and
% computes the estimation error between observed and true functional values.

% Args:
% sigma -- observer coefficient matrix: ord x (m + 2*p + ord)
% u_traj -- input trajectory: m x len matrix
% y_traj -- output trajectory: p x len matrix  
% z_traj -- true functional trajectory: r x len matrix

% Outputs:
% Z_error -- functional estimation error: r x len matrix

% Author: Ziyuan Luo
% email: ziyuan.luo@bit.edu.cn 
% Last revision: Oct-20-2025

%------------- BEGIN CODE --------------

    % Initialize observer parameters and state
    ord = size(sigma, 1);
    r = size(z_traj, 1);
    z0 = rand(ord, 1);
    len = size(u_traj, 2);
    z_all = z0;
    
    % Simulate functional observer over entire trajectory
    for t = 2:len
        z0 = sigma * [u_traj(:, t-1); y_traj(:, t-1); y_traj(:, t); z0];
        z_all = [z_all z0];
    end

    % Compute estimation error (first r components vs true functionals)
    Z_error = z_all(1:r, :) - z_traj;

%------------- END OF CODE --------------
end