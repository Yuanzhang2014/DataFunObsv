function X = SVD_subspace_intersection(A, B, tol1, tol2)
% SVD_SUBSPACE_INTERSECTION - Compute row space intersection using SVD
% Computes the intersection of row spaces of matrices A and B using
% singular value decomposition with economic computation for efficiency.

% Args:
% A -- first input matrix: row1 x T matrix
% B -- second input matrix: row2 x T matrix  
% tol1 -- (optional) tolerance for first SVD truncation
% tol2 -- (optional) tolerance for second SVD truncation

% Outputs:
% X -- matrix whose rows form a basis for the intersection of row spaces
%      of A and B: r x T matrix, where r = dimension of intersection

%------------- BEGIN CODE --------------

    % Form composite matrix and compute economic SVD
    H = [A; B];
    row = size(A, 1);
    [U, S, ~] = svd(H, 'econ');
    
    % Compute singular values and set first tolerance if not provided
    s1 = diag(S);
    if nargin == 2
        tol1 = max(size(H)) * eps(max(s1));
    end
    
    % Determine numerical rank and partition singular vectors
    r1 = sum(s1 > tol1);
    U11 = U(1:row, 1:r1);
    U12 = U(1:row, r1+1:end);
    S11 = S(1:r1, 1:r1);
    
    % Compute intermediate matrix and perform second economic SVD
    [UU, SS, ~] = svd(U12' * U11 * S11, 'econ');
    s2 = diag(SS);
    
    % Set second tolerance and extract intersection basis
    if nargin == 2
        tol2 = max(size(U12' * U11 * S11)) * eps(max(s2));
    end
    
    r2 = sum(s2 > tol2);
    Uq = UU(:, 1:r2);
    
    % Project back to original space for final intersection basis
    X = Uq' * U12' * A;

%------------- END OF CODE --------------
end
