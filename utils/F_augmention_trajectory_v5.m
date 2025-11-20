function [M1,zf,fullF,exist]=F_augmention_trajectory_v5(M1,zf,wp,wf,r) 
%  F_augmentation_trajectory_v5 -- augment functional state (find R ) 
%  such that rank([M1;R*xp]) == rank([M1;zf;R*xp;R*xf]) based on given 
%  trajectories xp xf

% Args: 
% M1 -- historical data matrix: [Up;Yp;Yf;Zp] 
% zf -- historical data matrix of future functional state
% wp & wf -- historical data matrix of available information
% r -- order of the original functional states (rank(L))

% Outputs:
% M1 -- augmented data matrix: [Up;Yp;Yf;Zp;fullF*wp]
% zf -- augmented data matrix of future functional state: [zf;fullF*wf]
% fullF -- augmented functional fullF = [L;R]
% exist -- indicator of whether all the conditions in the procedure 
% are satisfied

% Authors: Yuan Zhang, Ziyuan Luo
% email:  zhangyuan14@bit.edu.cn
%         ziyuan.luo@bit.edu.cn 
% Last revision: Oct-28-2025

%------------- BEGIN CODE --------------
exist=false;
n=size(wp,1);
fullF=zeros(0,n);
ite=0;
newF=[];

% find newF that increase rank(M1) but remain rank([M1;zf])
[M1,zf,newF]=F_augmention_v3(M1,zf,wp,wf);
while ( (check_consistent(M1,zf)==false) && rank_svd([M1;zf]) ~= rank(M1) ) && size(newF,1)<n-r
    M1=rowbasis(M1);
    T1=null_threshold([wp;M1]');
    T1=T1';
    Phip=T1(:,1:n);  % remains the rank_svd
    Phip_orth=null_threshold(Phip);  
    Phip_orth=Phip_orth'; % adding rank_svd one for the second matrix 
    if size(Phip_orth,1)>0
        tempsubspace=rowbasis([M1;zf]);
        T1=null_threshold([wp;tempsubspace]');
        T1=T1';
        Phip2=T1(:,1:n);        % remains the rank_svd
        T2=null_threshold([wf;tempsubspace]'); 
        T2=T2';
        Phif3=T2(:,1:n); % remains the rank_svd row space
        Phinewcase1=SVD_subspace_intersection(Phip2,Phip_orth)';
        if size(Phinewcase1,2)>0          % rk[M1;zf;f*xf] remains constant, rk[M1;zf;f*xp] increases.  
            M1=[M1;Phinewcase1'*wp];
            zf=[Phinewcase1'*wf;zf];
            newF=[newF;Phinewcase1'];
            [M1,zf,addF]=F_augmention_v3(M1,zf,wp,wf);
            newF=[newF;addF];
        else
            Phinewcase2=SVD_subspace_intersection(Phif3,Phip_orth)';
            if size(Phinewcase2,2)>0
                newF=[newF;Phinewcase2'];
                zf=[Phinewcase2'*wf;zf];
                M1=[M1;Phinewcase2'*wp]; 
                [M1,zf,addF]=F_augmention_v3(M1,zf,wp,wf);
                newF=[newF;addF];
            else
                n1=size(Phip_orth,1);
                n2=ceil(n1/3);
                newF=[newF;Phip_orth(1:n2,:)];
                zf=[Phip_orth(1:n2,:)*wf;zf];
                M1=[M1;Phip_orth(1:n2,:)*wp];
    
                %%  [M1,zf,newF]=F_augmention_v3(M1,zf,xp,xf);
                [M1,zf,addF]=F_augmention_v3(M1,zf,wp,wf);
                newF=[newF;addF];
            end      
       end 
    else
        % disp('Cannot find F saifying the equality');
        break; 
    end

    [M1,zf,addF]=F_augmention_v3(M1,zf,wp,wf);
    fullF=[fullF;newF;addF];
    fullF=rowbasis(fullF);
    % check whether the rank condition is satisfied
    if  check_consistent(M1,zf)
        exist=true;
        disp('A suitable F has been found');
        break;
    end
end
end


