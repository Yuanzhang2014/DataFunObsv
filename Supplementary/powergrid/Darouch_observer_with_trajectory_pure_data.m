function [Sigma, order] = Darouch_observer_with_trajectory_pure_data(u_traj, y_traj, z_traj)
% =========================================================================
% Filename: Darouch_observer_with_trajectory_pure_data.m
% Purpose: Design Darouch functional observer using pure trajectory data without system model
% Inputs:
%   u_traj - Input trajectory data (m × N)
%   y_traj - Output trajectory data (p1 × N)  
%   z_traj - Target functional state trajectory (r × N)
% Outputs:
%   Sigma - Observer parameter matrix
%   order - Boolean indicating if observer design was successful
% Main functions:
%   1. Construct data matrices from trajectory data
%   2. Perform singular value decomposition for numerical stability
%   3. Design functional observer using Darouch method
%   4. Check observability and pole placement conditions
% =========================================================================
Up = u_traj;
N = size(Up, 2);           
m = size(Up, 1);               
p1 = size(y_traj, 1);     
Sigma = [];
order = false;            

%% Data preparation section
up = Up(:, 1:end-1);          
yp = y_traj(:, 1:end-1);     
yf = y_traj(:, 2:end);          
zp = z_traj(:, 1:end-1);      
zf = z_traj(:, 2:end);        

p = [up; yp; yf; zp];

T = N - 1;
svd_p = svd(p * p' / T);   
svd_trun_index = find(svd_p > 1e-7);
[U1, ~, ~] = svd(p);  

Sigma = zf * pinv(p);
Sigmazp = Sigma(:, (m + p1 + p1 + 1):end);

%% Observer existence and stability conditions checking

if ((max(abs(eig(Sigmazp)))) <= 1 - 1e-1)
    order = true;
    disp('Order of the functional observer');
    disp(size(z_traj, 1));
    disp('type 1 - Naturally stable system');
    else
    nu = U1(:, max(svd_trun_index) + 1:end)';
    nuk = nu(:, (m + p1 + p1 + 1):end);
    desired_poles = 0.5 * (rand(size(Sigmazp, 1), 1) - 0.5);
        if size(nuk, 1) > 0 && (rank_svd(ctrb(Sigmazp', nuk')) == size(Sigmazp, 1))
            KK = place(Sigmazp', nuk', desired_poles);
            order = true;
            Sigma = Sigma - KK' * nu;
            disp('Order of the functional observer');
            disp(size(z_traj, 1));
            disp('type 2 - Stabilized by pole placement');
                else
                disp('The Darouch observer does not exist');
                Sigma = [];
                order = false;
                disp('type 3 - Observer does not exist');
        end
end