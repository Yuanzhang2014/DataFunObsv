% Learning feedback control laws using functional observers in load frequency control
% for LFC in the New England three area 39-bus Test case

%Data source:
%Shangguan, X. C., He, Y., Zhang, C. K., Jiang, L., & Wu, M. (2021). Adjustable event-triggered load frequency control of power systems using control-performance-standard-based fuzzy logic. IEEE Transactions on Fuzzy Systems, 30(8), 3297-3311.

                     
% Author: Wenxuan Xu
% email:  xuwenxuan1209@163.com
% Last revision: Oct-28-2025

%------------- BEGIN CODE --------------
clear;
close all;
addpath('../utils')
% Load grid parameters 
N = 3;                  % Number of control areas
n = 4;                  % State dimension of each area        

% Inertia constants of generators
M11=10; M12=6.06; M13=7.16;
M14=5.72; M15=5.20; M16=6.96; M17=5.28;
M18=4.86; M19=6.90; M20=8.40;

M1 = M11+M12+M13; 
M2 = M14+M15+M16+M17;  
M3 = M18+M19+M20;  

% Generator damping coefficien
D1 = 0.1; R1 = 0.05; Tt1 = 0.33; Tg1 = 0.1; beta1 = 2;
D2 = 0.1; R2 = 0.05; Tt2 = 0.25; Tg2 = 0.17; beta2 = 2;
D3 = 0.1; R3 = 0.05; Tt3 = 0.33; Tg3 = 0.2; beta3 = 2;

% Tie-line coupling coefficients
T12 = 0.4166; T13 = 1.3272; T23 = 0.2959;

w = 2 * pi;    % frequency

% Local linearized dynamic models for each area
A1 = [
    -D1/M1,   1/M1,     0,       -1/M1;
    0,      -1/Tt1,   1/Tt1,     0;
   -R1/Tg1,  0,      -1/Tg1,     0;
    0,       0,        0,        0
];

A2 = [
    -D2/M2,   1/M2,     0,       -1/M2;
    0,      -1/Tt2,   1/Tt2,     0;
   -R2/Tg2,  0,      -1/Tg2,     0;
    0,       0,        0,        0
];

A3 = [
    -D3/M3,   1/M3,     0,       -1/M3;
    0,      -1/Tt3,   1/Tt3,     0;
   -R3/Tg3,  0,      -1/Tg3,     0;
    0,       0,        0,        0
];

% Tie-line coupling matrices
A12 = zeros(4); A12(4,1) = -w * T12;
A13 = zeros(4); A13(4,1) = -w * T13;
A21 = zeros(4); A21(4,1) = -w * T12;
A23 = zeros(4); A23(4,1) = -w * T23;
A31 = zeros(4); A31(4,1) = -w * T13;
A32 = zeros(4); A32(4,1) = -w * T23;

% Full interconnected system matrix
A = [
    A1,  A12, A13;
    A21, A2,  A23;
    A31, A32, A3
];

% Control input matrices 
B1 = [0; 0; 1/Tg1; 0];
B2 = [0; 0; 1/Tg2; 0];
B3 = [0; 0; 1/Tg3; 0];
B = blkdiag(B1, B2, B3);

% Disturbance input matrices
F1 = [-1/M1; 0; 0; 0];
F2 = [-1/M2; 0; 0; 0];
F3 = [-1/M3; 0; 0; 0];
F = blkdiag(F1, F2, F3);

% Euler discretization of continuous-time system
dt = 0.01;          
n1 = size(A,1); 

Ad = eye(n1) + dt * A;  
Bd = dt * B;         
Dd = dt * F;           

% Output 
C1 = [beta1, 0, 0, 1];
C2 = [beta2, 0, 0, 1];
C3 = [beta3, 0, 0, 1];
C = blkdiag(C1, C2, C3);
Cd = dt * C;

% Augment system for integral action
Aa = [Ad, zeros(12, 3); 
      Cd, eye(3)];         

Ba = [Bd; zeros(3, 3)];    
Da = [Dd; zeros(3, 3)];    

Ca = [C, zeros(3, 3); 
      zeros(3,size(Cd,2)),eye(3)];

n2 = size(Aa,1);  
m = size(Ba,2);   % Number of control inputs

Cbar = Ca;

% Normalize Aa to stabilize numerical operations
Aa = Aa ./ (max(abs(eig(Aa))) * 1.0001);  

%% Generate simulation data for data-driven controller design
T = 5000;
X1 = zeros(n2,T);
UpALL=[]; XALL=[];
D = 0*ones(N,T);    % Zero external disturbance

% Random initial state
X1(:,1) = 0.1*[0.004;0.006;0.006;0.006;0.003;0.003;0.003;0.004;0.0071;0.006;0.009;0.005;0.006;0.009;0.005];

for k = 1:T
    Up(:,k) = -0.001*rand(m,1);       % Random control input
    X1(:,k+1) = Aa*X1(:,k) + Ba*Up(:,k) + Da*D(:,k);
end

UpALL = [UpALL,Up];
XALL = [XALL,X1(:,1:end-1)];

u_traj = UpALL;
x_traj1 = XALL;

%% Data-driven target output controller design
L1 = eye(n2);
z_traj = L1*x_traj1;

up = u_traj(:,1:end-1);
zp = z_traj(:,1:end-1);
zf = z_traj(:,2:end);

p0 = [up; zp; zf];
p1 = [up; zp];

% Check equivalence for controller existence
if rank(p0)==rank(p1)

    % Solve linear matrix inequality for controller design
    T = zf * pinv([up; zp]);
    T1=T(:,1:m); T2=T(:,m+1:end);

    setlmis([]); 
    P = lmivar(1, [size(T2,1), 1]);  
    K1 =lmivar(2, [size(T1,2), size(T2,1)]);  

    lmiterm([1 1 1 P], -0.99, 1);
    lmiterm([1 1 2 P], 1, T2');
    lmiterm([1 1 2 -K1], 1,T1');  
    lmiterm([1 2 2 P], -1,1);  

    lmis = getlmis;
    [tmin, xfeas] = feasp(lmis);

    P_sol = dec2mat(lmis, xfeas, P);
    K1_sol = dec2mat(lmis, xfeas, K1);

    K = K1_sol * inv(P_sol);
    L = K;
    r = size(L,1);
else 
    disp('Unable to design the target output controller')
end

%% Compute data-driven functional observer
[Sigma,flag,newL]=Reduced_order_observ_augmented_subspace_intersection_v3_partial(u_traj,x_traj1,x_traj1,Cbar,L,eye(size(L,2)));
%% Apply control for load frequency control (LFC)
T = 50000;     % End of region 3
T1 = 10000;    % End of region 1
T2 = 15000;    % End of region 2

X1 = zeros(n2,T+1);    % With control
X2 = zeros(n2,T+1);    % Without control

u_traj=[]; hu=[];

% Initial state for testing
X1(:,1) = 0*[0.004;0.006;0.006;0.006;0.003;0.003;0.003;0.004;0.0071;0.006;0.009;0.005;0.006;0.009;0.005];
X2(:,1) = X1(:,1);

r1 = size(newL,1);
u1 = 1*rand(r1,1);
u = u1(1:r,:);

D = 0*ones(N,T); 
D(1,T1+1:end) = 0.01;   % Strong disturbance applied to area 1

% Region 1: No control
for k = 1:T1  
    u = 0*rand(m,1);
    X1(:,k+1) = Aa*X1(:,k) + Ba*u + Da*D(:,k);
    X2(:,k+1) = Aa*X2(:,k) + Da*D(:,k);
    u_traj = [u_traj,u];
end

% Region 2: Still no feedback control
for k = T1+1:T2   
    u = 0*rand(m,1);
    X1(:,k+1) = Aa*X1(:,k)+Ba*u+Da*D(:,k);  
    X2(:,k+1) = Aa*X2(:,k)+Da*D(:,k);       
    u_traj = [u_traj,u];
end

% Region 3: Apply designed controller
for k =T2+1:T  
    X1(:,k+1) = Aa*X1(:,k)+Ba*u+Da*D(:,k);  
    X2(:,k+1) = Aa*X2(:,k)+Da*D(:,k);       

    Y1 = Cbar*X1(:,k+1);
    u1 = Sigma * [u; Cbar*X1(:,k); Y1; u1];  
    u = u1(1:r,:); 
    u_traj = [u_traj,u];
end

x_traj1 = X1(:,1:T);
x_traj2 = X2(:,1:T);

hu = L * x_traj1;      % Observer output
eu = u_traj - hu;      % Observer error

%% -----------------------------------------
%%       Plotting (All comments added)
%% -----------------------------------------

t = 0:dt:(T-1)*dt;  % Time vector

% Global font settings
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultTextFontSize', 10);
set(0, 'DefaultLegendFontSize', 10);

% Define region time boundaries
area1_end = T1*dt;
area2_end = T2*dt;
area3_end = max(t);

%% --- Figure 1: Area 1 frequency deviation ---
figure(1)

plot(0:dt:area1_end, x_traj1(1,1:T1+1), 'k', 'LineWidth',1); hold on
plot(area2_end:dt:max(t), x_traj1(1,T2+1:50000),'LineWidth',1,'Color', [44,162,95]/255); hold on
plot(area1_end:dt:max(t), x_traj2(1,T1+1:50000),'r','LineWidth',1);

% Add region background colors
yl = ylim; y_min = yl(1); y_max = yl(2);

fill([0, area1_end, area1_end, 0],...
      [y_min-0.005, y_min-0.005, y_max+0.015, y_max+0.015],...
      [230 230 230]/255,'EdgeColor','none','FaceAlpha',0.3);

fill([area1_end, area2_end, area2_end, area1_end],...
      [y_min-0.005, y_min-0.005, y_max+0.015, y_max+0.015],...
      [1 0.8 0.8],'EdgeColor','none','FaceAlpha',0.3);

fill([area2_end, area3_end, area3_end, area2_end],...
      [y_min-0.005, y_min-0.005, y_max+0.015, y_max+0.015],...
      [0.8 1 0.8],'EdgeColor','none','FaceAlpha',0.3);

% Mark event points
impact_idx = find(t >= 100, 1);
if ~isempty(impact_idx)
    plot(t(impact_idx), x_traj1(1, impact_idx), 'ko', 'MarkerSize', 8, 'LineWidth', 1.5);
end

control_idx = find(t >= 150, 1);
if ~isempty(control_idx)
    plot(t(control_idx), x_traj1(1, control_idx), 'ks', 'MarkerSize', 8, 'LineWidth', 1.5);
end

xlim([50,450]); ylim([-0.014,0.007]);
set(gca,'ytick',-0.01:0.005:0.005);
set(gca,'xtick',0:100:500);
set(gcf,'Position',[200 200 250 200]); 
ax=gca; ax.LineWidth=0.6; box on

%% --- Figure 2: Tie-line power of Area 1 ---
figure(2)
plot(0:dt:area1_end, x_traj1(4,1:T1+1), 'k','LineWidth',1); hold on
plot(area2_end:dt:max(t), x_traj1(4,T2+1:50000),'LineWidth',1,'Color',[44,162,95]/255); hold on
plot(area1_end:dt:max(t), x_traj2(4,T1+1:50000),'r','LineWidth',1);

% Background regions
yl = ylim; y_min=yl(1); y_max=yl(2);

fill([0, area1_end, area1_end, 0], [y_min-0.1 y_min-0.1 y_max+0.1 y_max+0.1],...
     [230 230 230]/255,'EdgeColor','none','FaceAlpha',0.3);

fill([area1_end area2_end area2_end area1_end], [y_min-0.1 y_min-0.1 y_max+0.1 y_max+0.1],...
     [1 0.8 0.8],'EdgeColor','none','FaceAlpha',0.3);

fill([area2_end area3_end area3_end area2_end], [y_min-0.1 y_min-0.1 y_max+0.1 y_max+0.1],...
     [0.8 1 0.8],'EdgeColor','none','FaceAlpha',0.3);

% Mark events
impact_idx2 = find(t >= 100, 1);
if ~isempty(impact_idx2)
    plot(t(impact_idx2), x_traj1(4, impact_idx2),'ko','MarkerSize',8,'LineWidth',1.5);
end

control_idx2 = find(t >= 150, 1);
if ~isempty(control_idx2)
    plot(t(control_idx2), x_traj1(4, control_idx2),'ks','MarkerSize',8,'LineWidth',1.5);
end

xlim([50,450]); set(gcf,'Position',[200 200 250 200]);
ax=gca; ax.LineWidth=0.6; box on

%% --- Figure 3: Control input & observer error ---
figure(3)

yyaxis left
plot(t, u_traj(1,:), 'LineWidth',1,'Color',[7,108,176]/255);
hold on;
plot(t, hu(1,:), '--','LineWidth',1,'Color',[7,108,176]/255);
ylabel('Control input / observer estimate');

set(gca,'YColor',[7,108,176]/255);
ylim([-0.6,3]);

yyaxis right
plot(t, eu(1,:), '--','LineWidth',1,'Color','r');
ylabel('Observer error');
set(gca,'YColor','r');

xlim([150 152]);
set(gca,'xtick',150.25:0.5:151.75);
ylim([-1,1]);
set(gca,'ytick',-0.8:0.4:0.8);

set(gcf,'Position',[200 200 250 200]); 
ax=gca; ax.LineWidth=0.6; box on

