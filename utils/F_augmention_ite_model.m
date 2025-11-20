function [T1, T2, T3, addF] = F_augmention_ite_model(T1, T2, T3, A)
% F_AUGMENTION_ITE_MODEL - Augment functional matrix via subspace intersection
% Performs one iteration of functional augmentation using model-based subspace
% intersection to expand the functional observer design.

% Args:
% T1 -- extended observability matrix [C; C*A]
% T2 -- combined matrix [T1; current functionals] 
% T3 -- extended matrix [T2; current functionals*A]
% A -- system matrix

% Outputs:
% T1, T2, T3 -- updated matrices with new functional augmentation
% addF -- newly added functional matrix

% Authors: Yuan Zhang, Ziyuan Luo
% email:  zhangyuan14@bit.edu.cn
%         ziyuan.luo@bit.edu.cn 
% Last revision: Oct-28-2025
%------------- BEGIN CODE --------------

    addF = [];
    
    % Check if augmentation is needed based on rank condition
    if rank_svd(T2) ~= rank_svd(T3)
        % Maintain row basis properties
        T1 = rowbasis(T1);
        T2 = rowbasis(T2);
        T3 = rowbasis(T3);
        
        % Compute null spaces for subspace intersection
        phi_psi = null_threshold([T3 * A; T3]')';
        n1 = size(T3, 1);
        phi = phi_psi(:, 1:n1);
        
        theta_xi = null_threshold([phi * T3; T2]')';
        n2 = size(phi, 1);
        theta = theta_xi(:, 1:n2);
        
        % Find orthogonal complement and compute intersection
        theta_orgh = null_threshold(theta * phi * T3);
        Gamma = SVD_subspace_intersection(theta_orgh', (phi * T3))';
        addF = Gamma';
        
        % Update matrices with new functional
        T2 = [T2; addF];
        T3 = [T3; addF; addF * A];
    end

%------------- END OF CODE --------------
end