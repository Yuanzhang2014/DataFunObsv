function [Phi_x, Phi_y] = lift_koopman(x, y, num_centers_x, num_centers_y, num_steps,centers_x,centers_y)
% =========================================================================
% Filename: lift_koopman.m
% Purpose: Perform Koopman operator lifting on state and output data using thin plate spline radial basis functions
% Inputs:
%   x - Original state trajectory data (state_dim × num_steps)
%   y - Original output trajectory data (output_dim × num_steps)
%   num_centers_x - Number of lifting dimensions for state data
%   num_centers_y - Number of lifting dimensions for output data
%   num_steps - Number of time steps in the trajectory
%   centers_x - (Optional) Pre-computed centers for state lifting (for testing phase)
%   centers_y - (Optional) Pre-computed centers for output lifting (for testing phase)
% Outputs:
%   Phi_x - Lifted state trajectory in Koopman space
%   Phi_y - Lifted output trajectory in Koopman space  
%   centers_x - Centers used for state radial basis functions
%   centers_y - Centers used for output radial basis functions
% Main functions:
%   Lift original state and output data to higher-dimensional space
% =========================================================================

    thin_plate_spline = @(x, c) max(norm(x-c), 1e-3)^2 * log(max(norm(x-c), 1e-3));
  %  centers_x = 1.5 * lhsdesign(num_centers_x, 2) - 0.75;
    basis_functions_x = cell(num_centers_x, 1);

    for i = 1:num_centers_x
        center0_x = centers_x(i, :);
        basis_functions_x{i} = @(x) thin_plate_spline(x, center0_x);
    end
    
  %  centers_y = 1.5 * lhsdesign(num_centers_y, 2) - 0.75;
    basis_functions = cell(num_centers_x, 1); 

    for i = 1:num_centers_y
        center0_y = centers_y(i, :);
        basis_functions_y{i} = @(x) thin_plate_spline(x, center0_y);  
    end
    
    phix = [];
    for i = 1:num_centers_x  
        for j = 1:num_steps 
            phix(i, j) = basis_functions_x{i}(x(:, j));
        end
    end
    
    phiy = [];
    for i = 1:num_centers_y 
        for j = 1:num_steps 
            phiy(i, j) = basis_functions_y{i}(y(:, j));
        end
    end
    
    Phi_x = [x; phix];
    Phi_y = [y; phiy];
    
    Phi_x = (Phi_x - mean(Phi_x, 2)) ./ std(Phi_x, 0, 2);
    Phi_y = (Phi_y - mean(Phi_y, 2)) ./ std(Phi_y, 0, 2);
end