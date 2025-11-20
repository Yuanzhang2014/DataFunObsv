function [sigma,eigdistri] = pure_IO_data_funobsvr_design_eigen_noise(u_traj,Y,Z,k,std)
%% This function designs reduced-order functional observers using noisy input-output data
%% std is the standard deviation of measurement noise

    

    T = size(u_traj,2);
    if nargin == 3
        k = 20; % t > l(A,C)  t must be greater than the observability index 
    end
    m = size(u_traj,1);
    p = size(Y,1);
    N = T - 2*k + 1;

    zp = Z(:,1:N-1);
    zf = Z(:,2:N);
    
    HankelU = Hankel(u_traj,2*k,N);
    Up = HankelU(1:m*k,:);
    % U1 = HankelU(1+m:m*(k+1),:);
    Uf = HankelU(m*k+1:2*m*k,:);
    % Uf = HankelU(m*(k+1)+1:end,:);
    
    HankelY = Hankel(Y,2*k,N);
    Yp = HankelY(1:p*k,:);
    % Y1 = HankelY(1+p:p*(k+1),:);
    Yf = HankelY(p*k+1:2*p*k,:);
    % Yf = HankelY(p*(k+1)+1:end,:);

    Wp = [Up;Yp];
    % W1 = [U1;Y1];
    Wf = [Uf;Yf];
    % Wf = [Uf;Yf];
    
    % Phip = subspace_intersection(W0',Wp',1e-4)';
    % Phif = subspace_intersection(W1',Wf',1e-4)';
    Phi = SVD_subspace_intersection(Wp,Wf);
    Phip = Phi(:,1:end-1);
    Phif = Phi(:,2:end);
   % disp('dimension of Phi');
    %disp(size(Phi));
    % Phif = SVD_subspace_intersection(W1,Wf);
    % Q = Phip(:,2:end) * pinv(Phif(:,1:end-1));
    % max(max(abs(Phip(:,2:end) - Q*Phif(:,1:end-1))))

    % up = u_traj(:,k+1:k+N);
    % yp = Y(:,k+1:k+N);
    % yf = Y(:,k+2:k+N+1);
    % zp = Z(:,k+1:k+N);
    % zf = Z(:,k+2:k+N+1);
    up = u_traj(:,k+1:k+N-1);
    yp = Y(:,k+1:k+N-1);
    yf = Y(:,k+2:k+N);
    zp = Z(:,k+1:k+N-1);
    zf = Z(:,k+2:k+N);

    M1 = [up;yp;yf;zp];
      T=N-1;
      
        svd_M1= svd(M1*M1'/T);
        eigdistri=svd_M1;
%         rownum=size(M1,1);
%         figure()
%         semilogy(1:rownum,svd_M1,'-*');
%           figure()
%         plot(1:rownum,svd_M1,'-*');
      
    if check_consistent(M1,zf) == false
        r=rank_svd(zp);

        % HpU = HankelU(:,1:N-1);
        % HfU = HankelU(:,2:N);
        % HpY = HankelY(:,1:N-1);
        % HfY = HankelY(:,2:N);
        % W = [HankelY*(eye(size(HankelU,2))-pinv(HankelU)*HankelU);HankelU];
            svd_M1= svd(M1*M1'/T);
   svd_trun_index=find(svd_M1>std^2*1.2);
   [U1,S1,V1] = svd(M1);
   M1trunc=U1(:,svd_trun_index)*S1(svd_trun_index,svd_trun_index)*V1(:,svd_trun_index)';
        
        [~,~,R,exist]=F_augmention_trajectory_v5(M1trunc,zf,Phip,Phif,r);
       % disp('R size');
      %  disp(size(R));
        
        
        
        M1 = [M1trunc;R*Phip];
        zf = [zf;R*Phif];
        disp('check rank condition')
        disp(check_consistent(M1,zf));
    end
    
    sigma = zf*pinv(M1);
    % sigma = zf\M1;
    
    sigmazp = sigma(:, end-size(sigma,1)+1:end);
    nu = null_threshold(M1',5e-5)';

    svd_M1= svd(M1*M1'/T);
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
        % desired_poles = 0.9 * rand(size(sigmazp, 1), 1);
        K = place(sigmazp', nu(:,end-size(sigma,1)+1:end)', desired_poles);
        % K = lmi_place(sigmazp',nu(:,end-size(sigma,1)+1:end)',0.9);
        sigma = sigma - K'*nu;
      %  disp('Order of the functional observer');
      %  disp(size(sigma,1));
     else
        disp('The Darouch observer does not exist');
        sigma = [];
    end

function Y=Hankel(X,t,N) % t is the number of block rows, N is the number of columns
    % Y=[];
    % for i=1:t
    %     Yi=[];
    %     for j=1:N
    %         Yi=[Yi,X(:,i+j-1)];
    %     end
    %     Y=[Y;Yi];
    % end
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