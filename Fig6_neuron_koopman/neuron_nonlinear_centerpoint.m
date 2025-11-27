% =========================================================================
% Filename: neuron_nonlinear.m
% Purpose: Generate nonlinear neuron network system trajectories with linearized comparison
% Inputs:
%   node_num - Number of neurons in the network
%   p - Output/sensor dimension
%   r - Dimension of states to be observed  
%   Topology - Network connection topology matrix
%   C - Output measurement matrix
%   L - Target state observation matrix
%   N - Number of time steps for simulation
%   X0 - Initial state vector
% Outputs:
%   X_non - Nonlinear system state trajectory
%   Y_non - Nonlinear system output trajectory
%   Z_non - Nonlinear system observed state trajectory
%   U - Input trajectory
% Main functions:
%   Generate input-output data pairs for system identification
% =========================================================================

function [X_non, Y_non, Z_non, U] = neuron_nonlinear(node_num, p, r, Topology, C, L, N, X0)

% Default: no noise (when fewer than 7 input parameters provided)
if nargin < 7
    sigma = 0; 
end

%% System dimension and parameter initialization
n = node_num * 3;               % System state dimension (3 states per neuron)
m = node_num;                   % Input dimension (each neuron has independent input)

% Ensure no self-connections in topology matrix
for i = 1:node_num
    Topology(i, i) = 0;
end

% Initialize system matrices for linearized model
A = zeros(n, n);                % State transition matrix
B = zeros(n, m);                % Input matrix  
X_affine = zeros(n, 1);         % Affine/constant term

% Input signal and trajectory initialization
U = 10 * rand(m, N);            % Random input sequence
X = zeros(n, N+1);              % Linear system state trajectory
X(:, 1) = X0;                   % Set initial condition
Y = zeros(p, N);                % Output trajectory
Z = zeros(r, N);                % Observed state trajectory

xop = rand(node_num, 1);        % Linearization operating point for Taylor expansion
dt = 0.0001;                    % Forward Euler discretization time step
X_affine0 = X_affine;
apm = 1;  

%% Nonlinear system simulation section
% Generate trajectory using full nonlinear neuron dynamics
X_non = zeros(n, N+1);
X_non(:, 1) = X(:, 1); 

for k = 1:N
    dx = zeros(n, 1);  
    
    for i = 1:node_num
        dx((i-1)*3+1, :) = -1 * X_non((i-1)*3+1, k)^3 + 3 * X_non((i-1)*3+1, k)^2 + ...
                           X_non((i-1)*3+2, k) - X_non((i-1)*3+3, k) + U(i, k);
        dx((i-1)*3+2, :) = 1 - 5 * X_non((i-1)*3+1, k)^2 - X_non((i-1)*3+2, k);
        dx((i-1)*3+3, :) = 0.001 * (X_non((i-1)*3+1, k) + 1.6 - X_non((i-1)*3+3, k));
        
        for j = 1:node_num
            if Topology(i, j) ~= 0
                dx((i-1)*3+1, :) = dx((i-1)*3+1, :) + X_non((j-1)*3+1, k) - X_non((i-1)*3+1);
            end
        end
    end

    X_non(:, k+1) = (X_non(:, k) + dt * dx) / apm;
end

Y_non = C * X_non(:, 1:N);      % Nonlinear output measurements
Z_non = L * X_non(:, 1:N);      % Nonlinear target states  
X_non = X_non(:, 1:N);          % Trim to N time steps for consistency

end
