function plot_95_envelope_v5(eigdistri_all, Z_error_all, N_test,std_val,SNR_all)
    % Plot 95% envelope for eigdistri and Z_error
    figure('Position', [100, 100, 1200, 800]);
    
    % ==================== eigdistri semilogy envelope plot ====================
    subplot(2, 1, 1);
    hold on;
    
    % Get maximum length of all eigdistri vectors
    max_length = 0;
    for i = 1:length(eigdistri_all)
        max_length = max(max_length, length(eigdistri_all{i}));
    end
    
    % Create matrix to store all data (pad with NaN for different lengths)
    eigdistri_matrix = nan(length(eigdistri_all), max_length);
    
    for i = 1:length(eigdistri_all)
        current_length = length(eigdistri_all{i});
        % Store original linear values (no log transformation)
        eigdistri_matrix(i, 1:current_length) = eigdistri_all{i}';
    end
    
    % Calculate 95% envelope using original linear values
    upper_envelope = zeros(1, max_length);
    lower_envelope = zeros(1, max_length);
    median_values = zeros(1, max_length);
    
    for j = 1:max_length
        valid_data = eigdistri_matrix(:, j);
        valid_data = valid_data(~isnan(valid_data)); % Remove NaN values
        
        if ~isempty(valid_data) && length(valid_data) > 5 % Ensure enough data points
            upper_envelope(j) = prctile(valid_data, 97.5);
            lower_envelope(j) = prctile(valid_data, 2.5);
            median_values(j) = median(valid_data);
        else
            upper_envelope(j) = NaN;
            lower_envelope(j) = NaN;
            median_values(j) = NaN;
        end
    end
    
    % Plot using semilogy with original linear data
    x_data = 1:max_length;
    valid_idx = ~isnan(upper_envelope) & ~isnan(lower_envelope);
    
    if any(valid_idx)
       x_valid = x_data(valid_idx);
    upper_valid = upper_envelope(valid_idx);
    lower_valid = lower_envelope(valid_idx);
      median_valid = median_values(valid_idx);

    % 1. 使用 plot 和 fill 绘制所有元素（此时纵轴仍是线性的）
    fill([x_valid, fliplr(x_valid)], [upper_valid, fliplr(lower_valid)], ...
         [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    plot(x_valid, upper_valid, 'b-+', 'LineWidth', 1.5);
    plot(x_valid, lower_valid, 'b-+', 'LineWidth', 1.5);
    plot(x_valid, median_valid, 'r-+', 'LineWidth', 2);

    for i = 1:min(4,length(eigdistri_all))
        current_data = eigdistri_matrix(i, valid_idx);
        plot(x_valid, current_data, 'Color', [0.5, 0.5, 0.5, 0.3], 'LineWidth', 0.5);
    end

    % 2. 在所有元素绘制完毕后，强制设置 Y 轴为对数刻度
    set(gca, 'YScale', 'log'); % 这是解决问题的关键行
    end
    % std_val,SNR_all
    xlabel('Eigenvalue index');
    ylabel('Eigenvalues');
    %title('95% Envelope of eigenvalues');
    snr_avg=mean(SNR_all);
    %
  %  title(sprintf('95%% Envelope of eigenvalues (\\sigma_{Y} = %.4f, SNR = %.4f)', std_val, snr_avg))
  multi_line_txt =sprintf('\\sigma_{Y} = %.2f, average SNR = %.2f', std_val, snr_avg);
    annotation('textbox', [0.7, 0.8, 0.2, 0.1], 'String', ...
             multi_line_txt, 'EdgeColor', 'none','FontSize', 11);
  %  legend('95% Envelope', '', '', 'Median', 'Individual Experiments', 'Location', 'best');
    grid on;
    
    % ==================== Z_error envelope plot ====================
    % Get Z_error dimensions
    r = size(Z_error_all{1}, 1);
    
    for row = 1:r
        subplot(2, r, r + row);
        hold on;
        
        % Extract all data for current row
        Z_row_data = zeros(length(Z_error_all), N_test);
        for i = 1:length(Z_error_all)
            Z_row_data(i, :) = Z_error_all{i}(row, :);
        end
        
        % Calculate 95% envelope
        upper_env = zeros(1, N_test);
        lower_env = zeros(1, N_test);
        median_vals = zeros(1, N_test);
        
        for t = 1:N_test
            upper_env(t) = prctile(Z_row_data(:, t), 97.5);
            lower_env(t) = prctile(Z_row_data(:, t), 2.5);
            median_vals(t) = median(Z_row_data(:, t));
        end
        
        % Plot envelope area
        time_points = 1:N_test;
        fill([time_points, fliplr(time_points)], [upper_env, fliplr(lower_env)], ...
             [1, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        plot(time_points, upper_env, 'r-', 'LineWidth', 1.5);
        plot(time_points, lower_env, 'r-', 'LineWidth', 1.5);
        
        % Plot median
     %   plot(time_points, median_vals, 'b-', 'LineWidth', 2);
        
        % Plot some representative experiment trajectories
        num_traces = min(4, length(Z_error_all));
        for i = 1:num_traces
            plot(time_points, Z_row_data(i, :), 'Color', [0.5, 0.5, 0.5, 0.2], 'LineWidth', 0.5);
              ylim([-2.5,10]);
        end
        
        xlabel('Time t');
      %  ylabel(['\e_', num2str(row), '(t)']);
        ylabel(['$e_{', num2str(row), '}(t)$'], 'Interpreter', 'latex');
     %   title(['95% Envelope of e(', num2str(row), ')']);
      %  legend('95% Envelope', '', '', 'Median', 'Experiment Traces', 'Location', 'best');
        grid on;
      
    end
    
    % Adjust subplot spacing
   % sgtitle('95% Envelope Analysis of 100 Experiments');
end