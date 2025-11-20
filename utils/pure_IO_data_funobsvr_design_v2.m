function sigma = pure_IO_data_funobsvr_design_v2(U,Y,Z,k)
% IO_data_funobsvr_design -- Utilize the proposed procedure to design a
% functional observer with Input, Output and Functional state data.

% Args: 
% U -- input trajectory: mxN matrix
% Y -- output trajectory: pxN matrix
% Z -- functional state trajectory rxN matrix
% k -- number of block rows, make sure k > observability index (system lag)
%      Default: k = round(n/p)+1

% Output: 
% sigma -- a stable coefficient matrix of the designed (possibly 
%          aumented) functional observer.                  
%          i.e., [z(t+1);r(t+1)] = sigma * [u(t);y(t);y(t+1);z(t);r(t)]

% Author: Ziyuan Luo
% email: ziyuan.luo@bit.edu.cn 
% Last revision: Oct-20-2025

%------------- BEGIN CODE --------------
    % Initialize
    % extract length of trajectory, dimensions of input & output and columns
    % of Hankel matrix
    T = size(U,2);
    if nargin == 3
        k = 20; % t > l(A,C)  t must be greater than the observability index 
    end
    m = size(U,1);
    p = size(Y,1);
    N = T - 2*k + 1; 
    
    % transform U,Y,Z trajactories into Hankel matrices form 
    HankelU = Hankel(U,2*k,N);
    Up = HankelU(1:m*k,:);
    Uf = HankelU(m*k+1:2*m*k,:);
    
    HankelY = Hankel(Y,2*k,N);
    Yp = HankelY(1:p*k,:);
    Yf = HankelY(p*k+1:2*p*k,:);

    Vp = [Up;Yp];
    Vf = [Uf;Yf];
    
    % compute Vp and Vf through subspace intersection based on SVD
    Phi = SVD_subspace_intersection(Vp,Vf);
    Phip = Phi(:,1:end-1);
    Phif = Phi(:,2:end);
    
    % extract corresponding input-output-functional state trajactories
    up = U(:,k+1:k+N-1);
    yp = Y(:,k+1:k+N-1);
    yf = Y(:,k+2:k+N);
    zp = Z(:,k+1:k+N-1);
    zf = Z(:,k+2:k+N);
    M1 = [up;yp;yf;zp];

    % functional observer design
    % check if condition (10) holds. If no, augment Zp and Zf with R*Phip 
    % and R*Phif first. If yes, compute sigma directly.
    if check_consistent(p,zf) == false
        r=rank_svd(zp);
        [~,~,R,exist]=F_augmention_trajectory_v5(M1,zf,Phip,Phif,r);
        
        M1 = [M1;R*Phip];
        zf = [zf;R*Phif];
    end
    
    % compute sigma
    sigma = zf*pinv(M1);
    
    % check if sigmazp is schur stable
    % if yes, return sigma
    % if no, check if corresponding null matrix nu is empty and pair 
    % (sigmazp, nu) is observable
    %   if yes, then apply pole placement to ensure stability of sigmazp
    %   if no, the functional observer of current order doesn't exist. Extra
    %   augmentation is needed.
    sigmazp = sigma(:, end-size(sigma,1)+1:end);
    svd_M1= svd(M1*M1'/(N-1));
    svd_trun_index=find(svd_M1>1e-8); % ER
    [U1,~,~] = svd(M1);
    nu=U1(:,max(svd_trun_index)+1:end)';      
    % nuk = nu(:, (m + p1 + p1 + 1):end);
    if nargin == 4
        desired_poles = 0.4*rand(size(sigma,1),1);
    end
    if max(abs(eig(sigmazp))) < 1 
        disp('Order of the functional observer');
        disp(size(sigma,1));
    elseif ~isempty(nu) && (rank_svd(ctrb(sigmazp', nu(:,end-size(sigma,1)+1:end)')) == size(sigmazp, 1))
        K = place(sigmazp', nu(:,end-size(sigma,1)+1:end)', desired_poles);
        sigma = sigma - K'*nu;
        disp('Order of the functional observer');
        disp(size(sigma,1));
     else
        disp('The Darouch observer does not exist');
        sigma = [];
    end

% constructin a Hankel matrix 
% Y = X_{0,t,N} -- t block rows and N columns
function Y=Hankel(X,t,N) 
    % t is the number of row blocks, N is the number of column blocks
    row = size(X(:,1),1);
    col = size(X(:,1),2);
    Y = zeros(t*row,N*col);
    for i = 1:t
        for j = 1:N
            Y((i-1)*row+1:i*row,(j-1)*col+1:j*col) = X(:,i+j-1);
        end
    end

end


end

%------------- END OF CODE --------------