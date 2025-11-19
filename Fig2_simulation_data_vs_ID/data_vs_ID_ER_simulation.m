%% SPARSE CASE WITH OBSERVABLE & UNOBSERVABLE SYSTEMS
% Comparison between data-driven functional observer and ID+Model-based 
% functional observer in ER (Erdos-Renyi) networks.

% Author: Ziyuan Luo
% email: ziyuan.luo@bit.edu.cn 
% Last revision: Oct-20-2025

clear
clc
addpath('../utils')

% Experimental parameters
n_lib = [50 60 70 80 90 100];      % System sizes to test
obsv_lib = [1 2];                  % 1: observable, 2: unobservable
degree = 10;                       % Network degree parameter
num_tries = 100;                   % Number of Monte Carlo simulations

% Initialize result arrays
timeFun1 = zeros(length(n_lib), num_tries, length(obsv_lib));  % Model-based timing
timeFun2 = zeros(length(n_lib), num_tries, length(obsv_lib));  % Data-driven timing
ord_m = zeros(length(n_lib), num_tries, length(obsv_lib));     % Model-based observer order
ord = zeros(length(n_lib), num_tries, length(obsv_lib));       % Data-driven observer order
RRMSE_m = zeros(length(n_lib), num_tries, length(obsv_lib));   % Model-based RRMSE
RRMSE = zeros(length(n_lib), num_tries, length(obsv_lib));     % Data-driven RRMSE

% Results file name with timestamp
str_save = ['comparision_ER_' num2str(num_tries) 'tries_' datestr(now, 'mm_dd_yy') '.mat'];
str_save = ['./' str_save];

%% Main experimental loop
for n_idx = 1:length(n_lib)
    for ob_case = obsv_lib
        for count = 1:num_tries
            % System parameters
            n = n_lib(n_idx);
            p_diag = 0.5;
            p_edge = degree / n;
            m = ceil(2 * n / 5);   % Number of inputs
            p = ceil(n / 2);       % Number of outputs  
            r = ceil(n / 10);      % Number of functionals
            len_data = ceil((2 * n + 1) * 3 / p) * (m + p + 1) + 100; % Minimum data length
            
            % Generate system with specified observability property
            while 1
                [A, B, C, L, ~, IsObsv] = sparse_ER_gen(n, m, p, r, p_edge, p_diag);
                if ob_case == 1 && IsObsv
                    break;  % Observable case
                end
                if ob_case == 2 && ~IsObsv  
                    break;  % Unobservable case
                end
            end
            
            % Generate training data
            [U_train, X_train] = model2traj(A, B, len_data, 'train');
            
            %% Model-based functional observer design
            tic
            [A_id, B_id, CL_id, ~, ~, ~] = sys_id_traj(U_train, X_train, [C * X_train; L * X_train]);
            if ~isempty(CL_id)
                C_id = CL_id(1:p, :);
                L_id = CL_id(p+1:p+r, :);
                [newL1, exist] = F_augmention_model_v1(A_id, C_id, L_id); 
            else
                newL1 = [];
            end
            
            if isempty(newL1)
                ord_m(n_idx, count, ob_case) = -1;
                timeFun1(n_idx, count, ob_case) = -1;
                Sigma1 = [];
            else
                disp('ID+model_based:')
                [Sigma1, order] = Darouch_observer_with_trajectory_v3(X_train, U_train, C_id, newL1);
                timeFun1(n_idx, count, ob_case) = toc;
                ord_m(n_idx, count, ob_case) = size(newL1, 1);
                disp('time cost:')
                disp(timeFun1(n_idx, count, ob_case))
            end
            
            %% Data-driven functional observer design  
            disp('Data-Driven:')
            tic
            Sigma2 = pure_IO_data_funobsvr_design_v2(U_train, C * X_train, L * X_train, 2 * round(n / p) + 1);
            ord(n_idx, count, ob_case) = size(Sigma2, 1);
            timeFun2(n_idx, count, ob_case) = toc;
            disp('time cost:')
            disp(timeFun2(n_idx, count, ob_case))
            
            %% Performance evaluation on test data
            [U_test, X_test] = model2traj(A, B, 2 * len_data, 'test');
            Z_test = L * X_test;
            
            % Model-based performance
            if isempty(Sigma1)
                RRMSE_m(n_idx, count, ob_case) = -1;
            else
                error_m = Simulation_agumentation_observer_with_trajectory_no_plot(Sigma1, X_test, U_test, C_id, newL1, L_id);
                acc_error_m = sqrt(trace(error_m(1:r, end-len_data+1:end)' * error_m(1:r, end-len_data+1:end)) / ...
                              trace(Z_test(1:r, end-len_data+1:end)' * Z_test(1:r, end-len_data+1:end)));
                RRMSE_m(n_idx, count, ob_case) = acc_error_m;
            end
            
            % Data-driven performance
            if isempty(Sigma2)
                RRMSE(n_idx, count, ob_case) = -1;
            else
                error_data = pure_data_simulation_no_plot(Sigma2, U_test, C * X_test, L * X_test);
                acc_error = sqrt(trace(error_data(1:r, end-len_data+1:end)' * error_data(1:r, end-len_data+1:end)) / ...
                             trace(Z_test(1:r, end-len_data+1:end)' * Z_test(1:r, end-len_data+1:end)));
                RRMSE(n_idx, count, ob_case) = acc_error;
            end
            
            % Periodic saving of results
            if mod(count, 10) == 0 || count == 1 
                save(str_save);
            end
        end
    end
end