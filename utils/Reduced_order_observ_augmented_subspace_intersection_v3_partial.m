function [Sigma,flag,newL]=Reduced_order_observ_augmented_subspace_intersection_v3_partial(U,X,W,C,L,W_transform,desired_poles)
%  Reduced_order_observ_augmented_subspace_intersection_v3_partial -- 
%  compute matrix Sigma with given (partial information) trajectories and augmented functional state 


% Args: 
% U -- historical input matrix 
% X -- historical state matrix 
% W -- historical partial state matrix 
% C & L -- matrices in state space model
% W_transform -- since partial state matrix W is a perturbation of X,
% augmented R has to postmultiply a transformation matrix
% desired_poles

% Outputs:
% Sigma -- coefficient matrix for Functional observer
% flag -- whether the observer exist
% newL -- augmented functional newL = [L;R]


% Authors: Yuan Zhang, Ziyuan Luo
% email:  zhangyuan14@bit.edu.cn
%         ziyuan.luo@bit.edu.cn 
% Last revision: Oct-30-2025

%------------- BEGIN CODE --------------

    n=size(X,1);
    % yes=data_functional_observability_test(U,C*X,L*X,n);
    Sigma=[];
    exist=false;
    newL=L;
    flag=0;
    Up=U;


    up = Up(:,1:end-1);
    xp = X(:, 1:end-1);
    xf = X(:, 2:end);
    wp = W(:,1:end-1);
    wf = W(:,2:end);
    yp = C * xp;
    yf = C * xf;
    zp = L * xp;
    zf = L * xf;
    %%%%%%%%%%%%%%%%%%% Observer design
    p = [up; yp; yf; zp];
    if check_consistent(p,zf)==false
        M1=p;
        % r=n-rank_svd(L);
        r = rank_svd(L);
        [~,~,addnewL,exist]=F_augmention_trajectory_v5(M1,zf,wp,wf,r);
        newL=[L;addnewL*W_transform];
         
    %%%%%%%%%%%%%%%%%%% Observer design
        if nargin == 7
            [Sigma,flag]=Darouch_observer_with_trajectory_v3(X,U,C,newL,desired_poles);
        else
            [Sigma,flag]=Darouch_observer_with_trajectory_v3(X,U,C,newL);
        end
            
        if exist
            exist=true;
            disp('Successful!')
        else
           disp('A suitable augmented has not been found');
        end 
    else
         [Sigma,flag]=Darouch_observer_with_trajectory_v3(X,U,C,newL);
          exist=true;
    end
%------------- END CODE --------------  
end


  


