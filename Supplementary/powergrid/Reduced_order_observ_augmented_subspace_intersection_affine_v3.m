function [Sigma,exist,newL]=Reduced_order_observ_augmented_subspace_intersection_affine_v3(u_traj,X,C,L)
n=size(X,1);
yes=data_functional_observability_test(u_traj,C*X,L*X,n);
% yes = 1;
Sigma=[];
exist=false;
newL=L;
flag=0;
if yes
    disp('This system is functionally observable from data');
else
       disp('Attension: This system is NOT functionally observable from data');
end
Up=u_traj;
if rank_svd([Up;X;ones(1,size(X,2))])==size([Up;X],1) + 1
    disp('The PE Condition is satisfied')
else
    disp('Attension: the PE Condition is not satisfied')
end
    up = Up(:,1:end-1);
    xp = X(:, 1:end-1);
    xf = X(:, 2:end);
    yp = C * xp;
    yf = C * xf;
    zp = L * xp;
    zf = L * xf;
    %%%%%%%%%%%%%%%%%%% Observer design
    p = [up; yp; yf; zp; ones(1,size(yp,2))];
    f = [p; zf];
    % if rank_svd(p)~=rank_svd(f)
         [Sigma,yes2]=Darouch_observer_with_trajectory_affine_v1(X,u_traj,C,newL);
         exist=true;
    end
%% Augmenting F by subspace intersection, M1=[Up',Yp',Yf',Zp']'    


  


