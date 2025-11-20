function [Sigma,order]=Darouch_observer_with_trajectory_pure_data(u_traj,y_traj,z_traj)
Up=u_traj;
N=size(Up,2);
m=size(Up,1); % number of inputs
p1 = size(y_traj, 1);
Sigma=[];
order= false;
    up = Up(:,1:end-1);
    yp = y_traj(:, 1:end-1);
    yf= y_traj(:, 2:end);
    zp = z_traj(:, 1:end-1);
    zf= z_traj(:, 2:end);
    %%%%%%%%%%%%%%%%%%% Observer design
     p = [up; yp; yf; zp];
     T=N-1;
     svd_p= svd(p*p'/T);
   svd_trun_index=find(svd_p>1e-7);
   [U1,~,~] = svd(p);
   Sigma=zf*pinv(p);
   Sigmazp = Sigma(:, (m + p1 + p1 + 1):end);
        nu=U1(:,max(svd_trun_index)+1:end)';
        nuk = nu(:, (m + p1 + p1 + 1):end);
        desired_poles = 0.5 * (rand(size(Sigmazp, 1), 1)-0.5);
       
   if ((max(abs(eig(Sigmazp)))) <= 1-1e-1)
        order= true;
        disp('Order of the functional observer');
        disp(size(L,1));
        disp('type 1');

   elseif size(nuk,1)>0 && (rank_svd(ctrb(Sigmazp', nuk')) == size(Sigmazp, 1))
        KK = place(Sigmazp', nuk', desired_poles);
        order= true;
        Sigma = Sigma - KK' * nu;
        
        disp('Order of the functional observer');
        disp(size(z_traj,1));
        disp('type 2');
   else
         disp('The Darouch observer does not exist');
         Sigma=[];
         order=false;
         disp('type 3');
   end

end