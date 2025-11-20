function [A,B,C,L,w,IsObsv] = sparse_BA_gen(n,m,p,r,m0,m_add,p_diag)

% sparse_BA_gen -- generates BA network with each node being a
% 2-dimensional subsystem

% Args:
% n -- number of nodes, then network dimension = 2n
% m -- dimension of inputs
% p -- dimension of outputs
% r -- dimension of functional states
% m0 -- parameter of BA network: initial nodes
% m_add -- parameter of BA network: each new node is connected to m_add 
%          nodes
% p_diag -- probability of self-loop, default p_diag=0.5

% Outputs:
% (A,B,C,L) -- state space model: x(t+1) = Ax(t) + Bu(t)
%                                 y(t) = Cx(t), z(t) = Lx(t)
% w -- adjacency matrix of the n-dimensional network (each node is a 2-dimensional system)
% IsObsv -- whether the generated system is observable: 
%           IsObsv=1 -> Observable
%           IsObsv=0 -> Unobservable

if nargin == 6
    p_diag = 0.5;
end

    
[w,~,~,~] = BA_network(n,m,p,r,m0,m_add);
w = ~~w .*rand(n);

A = zeros(2*n,2*n);
A(1:n,1:end) = [diag(0.5*rand(n,1)+0.5) 0.5*diag(rand(n,1)+0.5)].*[eye(n) eye(n)];
A(n+1:2*n,1:n) = diag(rand(n,1).*(rand(n,1) < p_diag));
A(n+1:2*n,n+1:2*n) = w + (diag(0.5*rand(n,1)+0.5)-sum(w,2)).*eye(n);
A = ~~A .* rand(2*n);
A = A / (1.01*max(abs(eig(A))));

bin = conncomp(digraph(A'));
if p >= max(bin)-1
    c1_idx = bin > 1;
    c_backup=find(bin(1:n) == 1);
    c_idx = c1_idx;
    c_idx(c_backup(randperm(length(c_backup),p-sum(c1_idx)))) = 1;
    IsObsv = 1; % observable
else
    c_idx = randperm(n,p);
    IsObsv = 0; % structral unobservable
end
idx = randperm(n,m+r);
b_idx = idx(:,1:m);
l_idx = idx(:,m+1:end);

% dedicate matrices
I = eye(n).*rand(n);
B = [zeros(n,m);I(:,b_idx)];
C = [I(c_idx,:) zeros(p,n)];
L = [zeros(r,n) I(l_idx,:)];
end