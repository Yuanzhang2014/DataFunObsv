function [U,X] = model2traj(A,B,len,mode)

% This function generate input and state sequences with state-space (A,B,C,L)
n = size(A,1);
m = size(B,2);

X = zeros(n,len);
if strcmp(mode, 'train') == 1
    U = 4*rand(m,len)-2;
else
    U = 2*sin(2*pi*rand(m,len));
    % U = ones(m,len);
end
X(:,1) = rand(n,1);

for k = 1:len-1
    X(:,k+1) = A*X(:,k) + B*U(:,k);
end

end