function [K, nu] = slack_place(Sigma, p, T)
% SLACK_PLACE - Stabilize observer via iterative pole placement
% Attempts to stabilize the functional observer by iteratively adjusting
% the null space selection for pole placement with desired pole constraints.

% Args:
% Sigma -- observer coefficient matrix
% p -- regression matrix used for null space computation
% T -- normalization factor for singular value calculation

% Outputs:
% K -- stabilization gain matrix
% nu -- selected null space basis used for stabilization

% Author: Ziyuan Luo
% Email: ziyuan.luo@bit.edu.cn
% Last revision: Nov.18, 2025
%------------- BEGIN CODE --------------

    % Initialize parameters and extract observer dynamics
    thresh = 1e-20;
    Sigmazp = Sigma(:, end-size(Sigma,1)+1:end);
    svd_p = svd(p * p' / T);
    
    % Generate desired poles: keep stable poles, replace unstable ones
    e = eig(Sigmazp);
    idx = abs(e) > 1;
    desired_poles = e;
    desired_poles(idx) = .8 * rand(sum(idx), 1) + 0.1;
    
    % Iterative pole placement with varying null space dimensions
    [U1, ~, ~] = svd(p);
    warning('off', 'all');
    
    for i = max(find(svd_p > thresh)) : -2 : 1
        nu = U1(:, i+1:end)';
        nuk = nu(:, end-size(Sigmazp,1)+1:end);
        
        % Attempt pole placement with current null space
        try 
            prec = [];
            [K, prec] = place(Sigmazp', nuk', desired_poles);
        catch ME
        end
        
        % Exit if successful pole placement achieved
        if prec > 0
            break
        end
    end
    
    warning('on', 'all');

%------------- END OF CODE --------------
end
