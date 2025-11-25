% This script generates an example for demo.m
% Attention: this code provides a customized procedure. For the example in
% out paper, please see 'demo.m' and load data file 'demo.mat'
% two System Identification models + Luenberger observer design with same
% poles
% data driven FO design with U-Y-Z data

% Author: Ziyuan Luo
% email:  ziyuan.luo@bit.edu.cn
% Last revision: Oct-28-2025

%------------- BEGIN CODE --------------
addpath('../utils') % Add the folder .\utils for functions
%% ID: 9-nodes system and 10-nodes system
n=10;
m=3;
p=2;
r=1;
p_edge = 0.08;
I = eye(n);
idx = randperm(n,m+p+r);
B = I(:,idx(1:m)) .* rand(n,m);
C = I(idx(m+1:m+p),:) .* rand(p,n);
L = I(idx(m+p+1:end),:) .* rand(r,n);
A = rand(n) < p_edge;
A = A + I(randperm(n,n),:);
A = A .* rand(n);
A = A / max(abs(eig(A))*1.1);

s = ceil ((n + 1) * 3 / p);
N = s * p + s * m;
T = N + s - 1;
X_train = zeros(n,T+1);
X_train(:,1) = rand(n,1);
U_train = 10*rand(m,T);
Y_train = zeros(p,T);
Z_train = zeros(r,T);
for k = 1 : T
    X_train(:,k+1) = A * X_train(:,k) + B * U_train(:,k);
    Y_train(:,k) = C * X_train(:,k);
    Z_train(:,k) = L * X_train(:,k);
end
[A_id1, B_id1, CL_id1, ~, ~, ~] = sys_id_traj(U_train, X_train, [Y_train;Z_train], 9);
C_id1 = CL_id1(1:p,:);
L_id1 = CL_id1(p+1:p+r,:);
[A_id2, B_id2, CL_id2, ~, ~, ~] = sys_id_traj(U_train, X_train, [Y_train;Z_train], 10);
C_id2= CL_id2(1:p,:);
L_id2 = CL_id2(p+1:p+r,:);

%% FO design
N=100;
% U_train = sin(rand(m,N));
U_train = 2*rand(m,N);
X_train = zeros(n,N);
X_train(:,1) = rand(n,1);
for i=1:N-1
    X_train(:,i+1) = A*X_train(:,i) + B*U_train(:,i);
end
Sigma = pure_IO_data_funobsvr_design_v2(U_train,C*X_train,L*X_train,n);


%% Luenberger Design

x0_id = zeros(n,1);
x0_id1 = zeros(size(A_id1,1),1);
x0_id2 = zeros(size(A_id2,1),1);

x0 = 5*rand(n,1);
x_real = x0;
x_id1 = x0_id1;
x_id2 = x0_id2;

z0 = zeros(size(Sigma,1),1);
z_all = z0;

desired_poles = eig(A);
% desired_poles = 0.6*rand(n,1);
K1 = place(A_id1',C_id1',[0.9091 + 0.0000i
  -0.4163 + 0.5780i
  -0.4163 - 0.5780i
  -0.0146 + 0.4096i
  -0.0146 - 0.4096i
   0.2780 + 0.0000i
  -0.2863 + 0.0000i
  -0.0309 + 0.0000i
  0.5642 + 0.0000i])'; % selected eigenvalues accroding to eig(A), example we use is provided in demo.mat
K2 = place(A_id2',C_id2',desired_poles)';
U_test = 3 + randn(m,N);

for i = 1:N-1
    x1 = A*x0 + B*U_test(:,i);
    y0 = C*x0;
    x1_id1 = A_id1*x0_id1 + B_id1*U_test(:,i) + K1*(y0 - C_id1 * x0_id1);
    x1_id2 = A_id2*x0_id2 + B_id2*U_test(:,i) + K2*(y0 - C_id2 * x0_id2);
    z1 = Sigma * [U_test(:,i);y0;C*x1;z0];
    
    x_real = [x_real x1];
    x_id1 = [x_id1 x1_id1];
    x_id2 = [x_id2 x1_id2];
    z_all = [z_all z1];
    
    x0 = x1;
    x0_id1 = x1_id1;
    x0_id2 = x1_id2;
    z0 = z1;
   
end


y_real = C*x_real;
z_real = L*x_real;
y_id1 = C_id1*x_id1;
z_id1 = L_id1*x_id1;
y_id2 = C_id2*x_id2;
z_id2 = L_id2*x_id2;

% save('./yourdemo.mat')

%------------- END OF CODE --------------
