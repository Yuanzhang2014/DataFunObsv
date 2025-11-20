function [M1,zf,fullF,exist]=F_augmention_trajectory_v5(M1,zf,xp,xf,r) 
%% Augmenting F by subspace intersection, M1=[Up',Yp',Yf',Zp']'
%% Here, r is the number of functional states M1,zf,xp,xf,r
exist=false;
%xp=rowbasis(xp);
n=size(xp,1);
fullF=zeros(0,n);
ite=0;
newF=[];
while exist==false && ite<20
    ite=ite+1;
   [M1,zf,newF]=F_augmention_v3(M1,zf,xp,xf);
%while rank_svd(M1)<rank_svd([M1;zf]) && size(newF,1)<n-r
while ( (check_consistent(M1,zf)==false) && rank([M1;zf]) ~= rank(M1) )&& size(newF,1)<n-r
     M1=rowbasis(M1);
     % T1=null_threshold([xp;M1]');
     T1=null([xp;M1]');
     T1=T1';
     Phip=T1(:,1:n);  % remains the rank_svd
     Phip_orth=null_threshold(Phip);  
     Phip_orth=Phip_orth'; % adding rank_svd one for the second matrix 
     if size(Phip_orth,1)>0
          tempsubspace=rowbasis([M1;zf]);
          % T1=null_threshold([xp;tempsubspace]');
          T1 = null([xp;tempsubspace]');
          T1=T1';
          Phip2=T1(:,1:n);        % remains the rank_svd
%           Phip_orth2=null_threshold(Phip2);
%           Phip_orth2=Phip_orth2'; % adding rank_svd one for the second matrix 
%           Phinew=SVD_subspace_intersection(Phip_orth2',Phip_orth');% column space
%       
           % T2=null_threshold([xf;tempsubspace]'); 
           T2=null([xf;tempsubspace]'); 
           T2=T2';
           Phif3=T2(:,1:n); % remains the rank_svd row space
%            Phif_orth3=null_threshold(Phif3); % adding rank_svd one for the third matrix 
           Phinewcase1=SVD_subspace_intersection(Phip2',Phip_orth');
           if size(Phinewcase1,2)>0          % rk[M1;zf;f*xf] remains constant, rk[M1;zf;f*xp] increases.  
               M1=[M1;Phinewcase1'*xp];
               zf=[Phinewcase1'*xf;zf];
               newF=[newF;Phinewcase1'];
             %  Zf=newF*Xf;
               [M1,zf,addF]=F_augmention_v3(M1,zf,xp,xf);
                newF=[newF;addF];
           else
               Phinewcase2=SVD_subspace_intersection(Phif3',Phip_orth');%  % rk[M1;zf;f*xp] remains constant, rk[M1;zf;f*xf] increases.
              % Phinewcase2=SVD_subspace_intersection(Phinew,Phif_orth3);
               if size(Phinewcase2,2)>0
                        newF=[newF;Phinewcase2'];
                        zf=[Phinewcase2'*xf;zf];
                        M1=[M1;Phinewcase2'*xp]; 
                        [M1,zf,addF]=F_augmention_v3(M1,zf,xp,xf);
                        newF=[newF;addF];
               else
%                     Phinewcase31=SVD_subspace_intersection(Phip_orth2',Phip_orth');% rk[M1;f*xp] increases, rk[M1;f*xf] increases.
%                     if size(Phinewcase31,2)>0
%                               newF=[newF;Phinewcase31'];
%                               zf=[Phinewcase31'*xf;zf];
%                           %   Zf=newF*Xf;
%                              M1=[M1;Phinewcase31'*xp];
%                   [M1,zf,addF]=F_augmention_v3(M1,zf,xp,xf);
%                   newF=[newF;addF];
%                     else       
%                         newF=[newF;Phip_orth(1,:)];
%                          zf=[Phip_orth(1,:)*xf;zf];
%                              M1=[M1;Phip_orth(1,:)*xp];
                             
                             %%%%%%%%%%%%%%
                             %%%%%%%%%%%%%
                             n1=size(Phip_orth,1);
                             n2=ceil(n1/3);
                                 newF=[newF;Phip_orth(1:n2,:)];
                         zf=[Phip_orth(1:n2,:)*xf;zf];
                             M1=[M1;Phip_orth(1:n2,:)*xp];
                              %%%%%%%%%%%%%%
                             %%%%%%%%%%%%%
                             %%  [M1,zf,newF]=F_augmention_v3(M1,zf,xp,xf);
                    [M1,zf,addF]=F_augmention_v3(M1,zf,xp,xf);
                   newF=[newF;addF];
               end      
           end 
     else
       %  disp('Cannot find F saifying the equality');
         break; 
     end
end
 
%if rank_svd(M1)<rank_svd([M1;zf]) && size(newF,1)==n-r
% if (check_consistent(M1,zf)==false) && size(newF,1)==n-r
% 
%   %   disp('Cannot find F saifying the equality');
% %else  
% end

[M1,zf,addF]=F_augmention_v3(M1,zf,xp,xf);
fullF=[fullF;newF;addF];
fullF=rowbasis(fullF);
%fullF=compact_firstrowL(fullF,L); 
%if rank_svd(M1)==rank_svd([M1;zf])
if  check_consistent(M1,zf)
    exist=true;
    disp('A suitable F has been found');
    break;
end
end
    end

