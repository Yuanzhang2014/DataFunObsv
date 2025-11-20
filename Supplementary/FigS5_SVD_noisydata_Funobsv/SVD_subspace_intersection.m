% function X = SVD_subspace_intersection(A,B,tol1,tol2)
% 
% This fuction realizes row space intersection for two matrices A and B 
% H = [A;B];
% row = size(A,1);
% 
% [U,S,~] = svd(H);
% s1 = S(S>0);
% if nargin == 2
%     tol1 = matlab.internal.math.getTolToCompareSVs(s1, max(size(H)));
% end
% r1 = sum(sum(s1>tol1));
% U11 = U(1:row,1:r1);
% U12 = U(1:row,r1+1:end);
% S11 = S(1:r1,1:r1);
% 
% [UU,SS,~] = svd(U12'*U11*S11);
% s2 = SS(SS>0);
% if nargin == 2
%     tol2 = matlab.internal.math.getTolToCompareSVs(s2, max(size(U12'*U11*S11)));
% end
% r2 = sum(sum(s2>tol2));
% Uq = UU(:,1:r2);
% 
% X = Uq'*U12'*A;
% 
% end
function X = SVD_subspace_intersection(A, B, tol1, tol2)

% This function realizes row space intersection for two matrices A and B
H = [A; B];
row = size(A, 1);

[U, S, ~] = svd(H, 'econ'); % 使用'econ'选项提高计算效率
s1 = diag(S); % 直接提取奇异值向量，更简洁
if nargin == 2
    % 替代方案：容差通常基于最大奇异值和矩阵维度设定
    tol1 = max(size(H)) * eps(max(s1));
end
r1 = sum(s1 > tol1); % 简化逻辑判断
U11 = U(1:row, 1:r1);
U12 = U(1:row, r1+1:end);
S11 = S(1:r1, 1:r1);

[UU, SS, ~] = svd(U12' * U11 * S11, 'econ');
s2 = diag(SS);
if nargin == 2
    % 同样的方法计算第二个容差
    tol2 = max(size(U12' * U11 * S11)) * eps(max(s2));
end
r2 = sum(s2 > tol2);
Uq = UU(:, 1:r2);

X = Uq' * U12' * A;
end