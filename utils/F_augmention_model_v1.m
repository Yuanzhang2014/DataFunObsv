function [newL, exist] = F_augmention_model_v1(A, C, L)
% F_AUGMENTION_MODEL_V1 - Functional augmentation for observer design
% Iteratively augments the functional matrix L to satisfy observer existence conditions
% using model-based approach with system matrices A, C, and initial functional L.

% Args:
% A -- system matrix: n x n
% C -- output matrix: p x n  
% L -- initial functional matrix: r x n

% Outputs:
% newL -- augmented functional matrix: (r + r_aug) x n
%         Contains original L as first rows plus augmented functionals
% exist -- boolean indicating if suitable augmentation was found
%          true: successful augmentation, false: failed to find augmentation

% Author: Yuan Zhang
% email: zhangyuan14@bit.edu.cn 
% Last revision: Oct-20-2025

%------------- BEGIN CODE --------------

    % Initialize variables
    exist = false;
    ite = 0;
    
    % Construct initial matrices for the iterative process
    T1 = [C; C * A];           % Extended observability matrix (2 steps)
    T2 = [T1; L];              % Combined observability and functionals
    T3 = [L * A; T2];          % One-step ahead functionals + T2
    
    newL = L;                  % Start with original functionals
    r = size(L, 1);
    n = size(A, 1);
    flag = 0;
    
    % Main augmentation loop (limited to 1 iteration)
    while exist == false && ite < 1
        ite = ite + 1;
        
        % Perform initial augmentation iteration
        [T1, T2, T3, addF] = F_augmention_ite_model(T1, T2, T3, A);
        newL = [newL; addF];
        
        % Iterative augmentation until conditions are satisfied
        count = 0;
        while ((check_consistent(T2, newL * A) == false) && ...
               rank_svd([T2; newL * A]) ~= rank(T2)) && ...
               size(newL, 1) < n - r && ...
               count <= 500
            
            % Maintain row basis properties
            T1 = rowbasis(T1);
            T2 = rowbasis(T2);
            T3 = rowbasis(T3);
            
            % Compute orthogonal complement for subspace intersection
            orth = null([T2; addF]);
            
            % Find new functionals through subspace intersection
            addF2 = SVD_subspace_intersection(T3, orth')';
            newL = [newL; addF2'];
            
            % Update matrices with new functionals
            T2 = [T2; addF2'];
            T3 = [T2; addF2' * A];
            
            % Perform another augmentation iteration
            [T1, T2, T3, addF] = F_augmention_ite_model(T1, T2, T3, A);
            
            % Add new functionals (skip if NaN values detected)
            if sum(sum(isnan(addF))) == 0
                newL = [newL; addF];
            end
            
            count = count + 1;
        end
        
        % Check for iteration limit failure
        if count >= 500
            newL = [];
            return;
        end
        
        % Final augmentation iteration
        [T1, T2, T3, addF] = F_augmention_ite_model(T1, T2, T3, A);
        newL = [newL; addF];
        
        % Check if augmentation conditions are satisfied
        if rank_svd([T2; newL * A]) == rank(T2)
            exist = true;
            disp('A suitable F has been found');
            break;
        end
    end
    
    % Ensure original L appears as first rows in the final basis
    newL = compact_firstrowL(newL, L);

%------------- END OF CODE --------------
end