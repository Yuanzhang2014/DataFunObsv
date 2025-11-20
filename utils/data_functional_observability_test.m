function yes = data_functional_observability_test(u_traj, y_traj, z_traj, n)
% DATA_FUNCTIONAL_OBSERVABILITY_TEST - Test functional observability from data
% Checks if the functional z can be determined from input u and output y data
% using a rank condition on Hankel matrices constructed from trajectories.

% Args:
% u_traj -- input trajectory: m x T matrix (m inputs, T time steps)
% y_traj -- output trajectory: p x T matrix (p outputs, T time steps)
% z_traj -- functional trajectory: r x T matrix (r functionals, T time steps)
% n -- system dimension (upper bound of observability index)

% Outputs:
% yes -- boolean indicating functional observability
%        true: functional is observable from input-output data
%        false: functional is not observable from input-output data

% Author: Ziyuan Luo
% email: ziyuan.luo@bit.edu.cn 
% Last revision: Oct-20-2025

%------------- BEGIN CODE --------------

    % Initialize result to false (not observable)
    yes = false;
    
    % Calculate trajectory parameters
    T = size(y_traj, 2);
    N = T + 1 - n;  % Number of columns in Hankel matrices
    
    % Construct Hankel matrices from input, output, and functional trajectories
    HankU = Hankel(u_traj, n, N);  % Input Hankel matrix
    HankY = Hankel(y_traj, n, N);  % Output Hankel matrix  
    HankZ = Hankel(z_traj, n, N);  % Functional Hankel matrix
    
    % Test functional observability using rank condition
    % The functional is observable if adding functional data doesn't increase
    % the rank beyond what's already captured by input and output data
    if rank([HankU; HankY; HankZ]) == rank([HankU; HankY])
        yes = true;
    end

%------------- END OF CODE --------------
end



