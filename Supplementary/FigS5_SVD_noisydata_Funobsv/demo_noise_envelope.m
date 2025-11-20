% Fixed parameter settings
% Generate Supplemental Figure S5: SVD-based dimension selection and functional observer performance with noisy measurement data.
n = 10;
m = 2;
p = 5;
r = 3;
N = 500;
N_test = 60;
num_experiments = 100; % Number of experiments

% Preallocate storage arrays
eigdistri_all = cell(1, num_experiments);
Z_error_all = cell(1, num_experiments);
SNR_all = zeros(1, num_experiments);

% Fixed part (generated outside the loop)
A = rand(n) < 0.6;
A = A .* rand(n);
A = A ./ (1.05 * max(abs(eig(A))));

I = eye(n);
idx = randperm(n, m + p + r);
b_idx = idx(1:m);
c_idx = idx(m+1:m+p);
l_idx = idx(m+p+1:end);

B = I(:, b_idx);
C = I(c_idx, :);
L = I(l_idx, :);

% Run 100 experiments
% load('noisesystem2.mat');
for exp_idx = 1:num_experiments
    fprintf('Running experiment %d/%d\n', exp_idx, num_experiments);
    
    % Random part (regenerated for each experiment)
    X_train = zeros(n, N);
    U_train = rand(m, N); % Uniform distribution
    X_train(:, 1) = rand(n, 1);

    for k = 1:N-1
        X_train(:, k+1) = A * X_train(:, k) + B * U_train(:, k);
    end
    
    Y_train = C * X_train;
    
    %% Add measurement noise
    std_val = 1; % Change the standard deviation to be 0, 0.01, 0.1, 1
    yfnoise = normrnd(0, std_val, size(Y_train));
    SNR = 10*log10(norm(Y_train, 'fro')^2 / (norm(yfnoise, 'fro')^2 + 10e-6));
    SNR_all(exp_idx) = SNR;
    Y_train = yfnoise + Y_train;

    %% Call subroutine
    [sigma, eigdistri] = pure_IO_data_funobsvr_design_eigen_noise(...
        U_train, Y_train, L * X_train, round(n/p) + 1, std_val);

    % Test data generation
    X_test = zeros(n, N_test);
    U_test = 10 * rand(m, N_test); % Random input
    X_test(:, 1) = (rand(n, 1) - 0.5);

    for k = 1:N_test-1
        X_test(:, k+1) = A * X_test(:, k) + B * U_test(:, k);
    end

    Y_test = C * X_test;
    Z_test = L * X_test;

    % Simulation
    ord = size(sigma, 1);
    r = size(Z_test, 1);
    z0 = 10 * rand(ord, 1);
    % z0 = zeros(ord, 1);
    z_all = z0;
    
    for t = 2:N_test
        z0 = sigma * [U_test(:, t-1); Y_test(:, t-1); Y_test(:, t); z0];
        z_all = [z_all, z0];
    end

    Z_error = z_all(1:r, :) - Z_test;
    
    % Store results
    eigdistri_all{exp_idx} = eigdistri;
    Z_error_all{exp_idx} = Z_error;
end

% Plot 95% envelope
plot_95_envelope_v5(eigdistri_all, Z_error_all, N_test, std_val, SNR_all);
% plot_95_envelope_v6(eigdistri_all, Z_error_all, N_test, std_val, SNR_all);