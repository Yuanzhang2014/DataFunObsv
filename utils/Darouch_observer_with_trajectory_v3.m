function [Sigma, flag] = Darouch_observer_with_trajectory_v3(X, U, C, L, desired_poles)
% DAROUCH_OBSERVER_WITH_TRAJECTORY_V3 - Design Darouch functional observer from trajectory data
% Designs a functional observer for the system using input-state-output trajectories

% Args:
% X -- state trajectory: n x N matrix (n states, N time steps)
% U -- input trajectory: m x N matrix (m inputs, N time steps)
% C -- output matrix: p1 x n matrix (p1 outputs, n states)
% L -- functional matrix: r x n matrix (r functionals to estimate)
% desired_poles -- (optional) desired poles for observer stabilization

% Outputs:
% Sigma -- observer coefficient matrix: r x (m + 2*p1 + r) matrix
%          Observer equation: z_hat(t+1) = Sigma * [u(t); y(t); y(t+1); z_hat(t)]
% flag -- status indicator:
%         1: observer naturally stable
%         0: observer stabilized with slack_place method
%        -1: slack_place method had large error
%        -2: insufficient null space for stabilization

% Author: Ziyuan Luo
% email: ziyuan.luo@bit.edu.cn 
% Last revision: Oct-20-2025

% Initialize outputs
Sigma = [];
flag = false;

% Check persistence of excitation condition
Up = U;
N = size(Up, 2);
m = size(Up, 1);
p1 = size(C, 1);


if rank_svd([Up; X]) == size([Up; X], 1)
    % PE condition satisfied - proceed with observer design
else
    disp('Attention: the PE Condition is not satisfied')
end

% Extract time-shifted trajectories for observer design
up = Up(:, 1:end-1);
xp = X(:, 1:end-1);
xf = X(:, 2:end);
yp = C * xp;
yf = C * xf;
zp = L * xp;
zf = L * xf;

% Construct regression matrix and compute observer coefficients
p = [up; yp; yf; zp];
Sigma = zf * pinv(p);
Sigmazp = Sigma(:, (m + p1 + p1 + 1):end);

% Observer stability analysis and stabilization
if nargin == 4
    % Check if observer is naturally stable
    if ((max(abs(eig(Sigmazp)))) < 1)
        disp('Order of the functional observer');
        disp(size(L, 1));
        flag = 1;
    else
        % Observer unstable - attempt stabilization using null space
        T = N - 1;
        svd_p = svd(p * p' / T);
        svd_trun_index = find(svd_p > 5e-9);
        [U1, ~, ~] = svd(p);
        nu = U1(:, max(svd_trun_index) + 1:end)';
        nuk = nu(:, (m + p1 + p1 + 1):end);
        
        if size(nuk, 1) > 0
            % Use slack_place method for stabilization
            [KK, nu] = slack_place(Sigma, p, T);
            Sigma = Sigma - KK' * nu;
            error = max(max(abs(KK' * nu * p)));
            
            if error < 0.07
                flag = 0;  % Successful stabilization
            else
                flag = -1; % Large stabilization error
            end
        else 
            flag = -2; % Insufficient null space
        end
        disp('Order of the functional observer');
        disp(size(L, 1));
    end
end

% Pole placement stabilization when desired poles are provided
if nargin == 5
    if size(nuk, 1) > 0 && ((rank_svd(ctrb(Sigmazp', nuk')) == size(Sigmazp, 1)))
        KK = place(Sigmazp', nuk', desired_poles(1:size(Sigmazp, 1)));
        Sigma = Sigma - KK' * nu;
        disp('Order of the functional observer');
        disp(size(L, 1));
    else
        disp('The Darouch observer does not exist');
        Sigma = [];
        flag = false;
    end
end

end