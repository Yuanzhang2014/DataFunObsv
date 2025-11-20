function [A,B,C,L,w,IsObsv] = sparse_ER_gen(n,m,p,r,p_threshold,p_diag)

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

% Author: Ziyuan Luo
% Email: ziyuan.luo@bit.edu.cn
% Last revision: Nov.18, 2025

if nargin == 5
    p_diag = 0.6;
end

    
w = rand(n) < p_threshold;
w = w.*~(eye(n)).*rand(n);
% w = w - w.*eye(n) + diag(rand(n,1).*(rand(n,1) < p_diag));
A = zeros(2*n,2*n);
A(1:n,1:end) = [diag(0.5*rand(n,1)+0.5) 0.5*diag(rand(n,1)+0.5)].*[eye(n) eye(n)];
A(n+1:2*n,1:n) = diag(rand(n,1).*(rand(n,1) < p_diag));
A(n+1:2*n,n+1:2*n) = w;
A = ~~A .* rand(2*n);
A = A / (1.01*max(abs(eig(A))));

c_idx = zeros(1,2*n);
[bin, binsize] = conncomp(digraph(A'));
[main_size, main_class] = max(binsize);
if p >= 2*n - main_size
    for i=1:max(bin)
        if i~= main_class
            tmp_idx = find(bin(1:n) ==i);
            c_idx(tmp_idx(1)) = 1;
        end
    end
    main_class_idx = find(bin(1:n) == main_class);
    c_idx(main_class_idx(randperm(length(main_class_idx),p-sum(c_idx)))) = 1;
    IsObsv = 1; % observable
else
    c_idx = randperm(n,p);
    IsObsv = 0; % structral unobservable
end
idx = randperm(n,m+r);
b_idx = idx(:,1:m);
l_idx = idx(:,m+1:end);


I = eye(n).*rand(n);
B = [zeros(n,m);I(:,b_idx)];
C = [I(~~c_idx,:) zeros(p,n)];
L = [zeros(r,n) I(l_idx,:)];

% u_traj = 4*rand(m,len_data)-2;
% X = zeros(2*n,len_data);
% X(:,1) = rand(2*n,1);
% for i=1:len_data-1
%     X(:,i+1) = A * X(:,i) + B * u_traj(:,i);
% end

% rank(ctrb(A,B))
% rank(obsv(A,C))
% rank([obsv(A,C);L])





end