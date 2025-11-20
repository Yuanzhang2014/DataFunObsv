function Z_error = Simulation_agumentation_observer_with_trajectory_no_plot(Sigma, x_traj, u_traj, C, L, original, mode)
% SIMULATION_AGUMENTATION_OBSERVER_WITH_TRAJECTORY_NO_PLOT - Simulate augmented observer performance
% Simulates the augmented functional observer using the coefficient matrix Sigma
% and computes the estimation error between observed and true functional values.

% Args:
% Sigma -- augmented observer coefficient matrix
% x_traj -- state trajectory: n x len matrix
% u_traj -- input trajectory: m x len matrix
% C -- output matrix: p x n
% L -- functional matrix: r x n
% original -- original functional matrix (unused in current implementation)
% mode -- simulation mode: 'acc' for accurate initialization, 'stochastic' for random

% Outputs:
% Z_error -- functional estimation error: r x len matrix

% Author: Ziyuan Luo
% email: ziyuan.luo@bit.edu.cn 
% Last revision: Oct-20-2025

%------------- BEGIN CODE --------------

    % Initialize simulation parameters and modes
    if nargin == 6
        mode = 'stochastic';
    end
    
    Len_data = size(x_traj, 2);
    X = x_traj;
    X0 = X(:, 1);
    Y0 = C * X0;
    r = size(Sigma, 1);
    
    % Set observer initial condition based on mode
    if strcmp(mode, 'acc') == 1
        observer_Z0 = L * X0;
    else
        observer_Z0 = rand(r, 1);
    end
    
    observer_Zall = observer_Z0;
    Zall = L * X;
    
    % Simulate augmented observer over trajectory
    for t = 2:size(u_traj, 2)
        U = u_traj(:, t-1);
        Y1 = C * X(:, t);
        Z1 = L * X(:, t);
        
        % Observer update using Sigma matrix
        observer_Z1 = Sigma * [U; Y0; Y1; observer_Z0];
        observer_Z0 = observer_Z1;
        Y0 = Y1;
        observer_Zall = [observer_Zall observer_Z1];
    end
    
    % Compute estimation error between true and observed functionals
    Z_error = Zall - observer_Zall;

%------------- END OF CODE --------------
end