function X = SVD_subspace_intersection(A, B, tol1, tol2)
% SVD_SUBSPACE_INTERSECTION - Compute row space intersection using SVD
% Computes the intersection of row spaces of matrices A and B using
% singular value decomposition for numerical stability and accuracy.

% Args:
% A -- first input matrix: row1 x T matrix
% B -- second input matrix: row2 x T matrix  
% tol1 -- (optional) tolerance for first SVD truncation
% tol2 -- (optional) tolerance for second SVD truncation

% Outputs:
% X -- matrix whose rows form a basis for the intersection of row spaces
%      of A and B: r x T matrix, where r = dimension of intersection

% Reference: Moonen et al.  On- and off-line identification of linear 
% state-space models[J]. International journal of control, 1989, 49(1): 219-232.

% Author: Ziyuan Luo
% Email: ziyuan.luo@bit.edu.cn
% Last revision: Oct.20, 2025

%------------- BEGIN CODE --------------

    % Form composite matrix and compute initial SVD decomposition
    H = [A; B];
    row = size(A, 1);
    T = size(A, 2) - 1;
    [U, S, ~] = svd(H);
    
    % Compute singular values and set first tolerance if not provided
    s1 = svd(1/(T-1) * (H * H'));
    if nargin == 2
        tol1 = matlab.internal.math.getTolToCompareSVs(s1, max(size(H)));
    end
    
    % Determine numerical rank and partition singular vectors
    r1 = sum(sum(s1 > tol1));
    U11 = U(1:row, 1:r1);
    U12 = U(1:row, r1+1:end);
    S11 = S(1:r1, 1:r1);
    
    % Compute intermediate matrix and perform second SVD
    H1 = U12' * U11 * S11;
    [UU, SS, ~] = svd(H1);
    
    % Set second tolerance and extract intersection basis
    s2 = svd(1/(T-1) * (H1 * H1'));
    if nargin <= 3
        tol2 = matlab.internal.math.getTolToCompareSVs(s2, max(size(U12' * U11 * S11)));
    end
    
    r2 = sum(sum(s2 > tol2));
    Uq = UU(:, 1:r2);
    
    % Project back to original space for final intersection basis
    X = Uq' * U12' * A;

%------------- END OF CODE --------------
end