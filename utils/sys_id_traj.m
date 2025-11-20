function [A_id, B_id, C_id, x0_id, Y, Y_id] = sys_id_traj(U, X, Y, n)
% SYS_ID_TRAJ - System identification from input-state-output trajectories
% Identifies system matrices A, B, C and initial state x0 using subspace methods

% Args:
% U -- input trajectory: m x T matrix (m inputs, T time steps)
% X -- state trajectory: nn x (T+1) matrix (nn states, T+1 time steps)
% Y -- output trajectory: p x T matrix (p outputs, T time steps)  
% n -- (optional) system order (number of states). If not provided, uses size(X,1)

% Outputs:
% A_id -- identified state transition matrix: n x n
% B_id -- identified input matrix: n x m
% C_id -- identified output matrix: p x n
% x0_id -- identified initial state: n x 1
% Y -- original output trajectory (passed through)
% Y_id -- simulated output from identified system: p x T

% Author: Ziyuan Luo
% email: ziyuan.luo@bit.edu.cn 
% Last revision: Oct-20-2025

% Initialize outputs
A_id = [];
B_id = [];
C_id = [];
x0_id = [];
Y_id = [];

% System dimensions and parameter setup
nn = size(X, 1);
m = size(U, 1);
p = size(Y, 1);
s = ceil((nn + 1) * 3 / p);  % Number of block rows for Hankel matrix
N = s * p + s * m;           % Number of columns for Hankel matrix  
T = N + s - 1;               % Required trajectory length

% Check data length sufficiency
if size(X, 2) < T
    disp('Warning: Insufficient Length of Data')
    return
end

% Handle optional system order parameter
if nargin == 4
    nn = n;
end

% Construct Hankel matrices for subspace identification
HankU = Hankel(U, s, N);
HankY = Hankel(Y, s, N);

% Check for numerical issues
if sum(isnan(HankY), 'all') ~= 0 || sum(isinf(HankY), 'all') ~= 0
    disp('Numerical ERROR -- NAN or INF')
    return
end

disp('sys_id is PE');

% Subspace identification using QR decomposition and SVD
[~, R] = qr([HankU; HankY]');
R = R';
R22 = R(s*m+1:end, s*m+1:s*m+s*p);
[Un, ~, ~] = svd(R22);

% Extract identified C and A matrices
C_id = Un(1:p, 1:nn);
J = Un(1:(s-1)*p, 1:nn);
K = Un(p+1:end, 1:nn);
A_id = pinv(J) * K;

% Least squares identification of initial state and B matrix
len = ceil((nn * m * 3) / p);
Y_col = zeros(len * p, 1);
phi = zeros(p * len, nn + nn * m);

for k = 1:len
    % Build output vector and observability matrix
    Y_col((k-1)*p+1:k*p, 1) = Y(:, k);
    if k == 1
        phi((k-1)*p+1:k*p, 1:nn) = C_id;
    else
        phi((k-1)*p+1:k*p, 1:nn) = phi((k-2)*p+1:(k-1)*p, 1:nn) * A_id;
    end
    
    % Build controllability/input influence matrix
    H = zeros(p, nn * m);
    for i = k-1:-1:1
        if i == k-1
            M = C_id;
        else
            M = M * A_id;
        end
        L = zeros(p, m * nn);
        for j = 1:m
            L(:, (j-1)*nn+1:j*nn) = U(j, i)' * M;
        end
        H = H + L;
    end
    phi((k-1)*p+1:k*p, nn+1:end) = H;
end

% Solve least squares problem for initial state and B matrix
Theta = pinv(phi) * Y_col;
x0_id = Theta(1:nn, :);
B_id = zeros(nn, m);

for i = 1:m
    B_id(:, i) = Theta((i-1)*nn+1+nn:i*nn+nn, :);
end

% Simulate identified system
Y_id = zeros(p, T);
X_id = zeros(nn, T+1);
X_id(:, 1) = x0_id;

for k = 1:T
    X_id(:, k+1) = A_id * X_id(:, k) + B_id * U(:, k);
    Y_id(:, k) = C_id * X_id(:, k);
end

end

function Y = Hankel(X, t, N)
% HANKEL - Construct Hankel matrix from trajectory data
% Args: X - input trajectory matrix, t - block rows, N - columns
    p = size(X, 1);
    Y = zeros(t * p, N);
    for i = 1:N
        for j = 1:t
            Y((j-1)*p+1:j*p, i) = X(:, i+j-1);
        end
    end
end