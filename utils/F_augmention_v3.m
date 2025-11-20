function [M1, zf, newF] = F_augmention_v3(M1, zf, xp, xf)
% F_AUGMENTION_V3 - Augment functional matrix via data-driven subspace intersection
% Performs one iteration of functional augmentation using trajectory data to
% expand the functional observer design while maintaining consistency conditions.

% Args:
% M1 -- regression matrix [Up'; Yp'; Yf'; Zp']
% zf -- future functional values
% xp -- past state trajectory
% xf -- future state trajectory

% Outputs:
% M1 -- updated regression matrix with new functionals
% zf -- updated future functional values
% newF -- newly added functional matrix

% Author: Yuan Zhang
% email: zhangyuan14@bit.edu.cn
% Last revision: Oct-20-2025

%------------- BEGIN CODE --------------

    newF = [];
    M1 = rowbasis(M1);
    
    % Check if augmentation is needed based on consistency condition
    if check_consistent(M1, zf) == false
        n = size(xp, 1);
        
        % Compute null spaces for state and data matrices
        T1 = null_threshold([xp; rowbasis([M1; zf])]')';
        Phip = T1(:, 1:n);
        
        T2 = null_threshold([xf; rowbasis([M1; zf])]')';
        Phif = T2(:, 1:n);
        
        % Find intersection subspace between past and future
        Phir = SVD_subspace_intersection(Phip, Phif)';
        Phir = Phir';
        
        % Compute orthogonal complement for augmentation
        T3 = null_threshold([xp; rowbasis(M1)]')';
        Phil = T3(:, 1:n);
        
        nu = null(Phil);
        
        % Find new functionals through subspace intersection
        Fadd = SVD_subspace_intersection(nu', Phir)';
        
        % Update matrices with new functionals
        newF = [newF; Fadd'];
        M1 = [M1; newF * xp];
        M1 = rowbasis(M1);
        zf = rowbasis([zf; Fadd' * xf]);
    end

%------------- END OF CODE --------------
end