function [Sigma,order]=Darouch_observer_with_trajectory_affine(x_traj,u_traj,C,L)
Up=u_traj;
N=size(Up,2);
m=size(Up,1); % number of inputs
p1 = size(C, 1);
x=x_traj;
Sigma=[];
order= false;
if rank_svd([Up;x;ones(1,size(Up,2))])==size([Up;x],1)+1
    disp('The PE Condition is satisfied')
else
    disp('Attention: the PE Condition is not satisfied')
end

    up = Up(:,1:end-1);
    xp = x(:, 1:end-1);
    xf = x(:, 2:end);
    yp = C * xp;
    yf = C * xf;
 
    %% adding measurement noises
  %%  yfnoise=normrnd(0,sigma,size(yp));
 %%   yp=yp+yfnoise;
 %%   yf(:,1:end-1)=yp(:,2:end);
 %%   yf(:,end)=yf(:,end)+normrnd(0,sigma,size(yp,1),1);
    zp = L * xp;
    zf = L * xf;
    %%%%%%%%%%%%%%%%%%% Observer design
    p = [up; yp; yf; zp; ones(1,size(yp,2))];
   % f = [p; zf];
%if rank_svd(p)==rank_svd(f)   
 % if check_consistent(p,zf)   
 if 1
   % Sigma = zf * (pinv(p));
   %%% corrected by noise
     T=N-1;
     svd_p= svd(p*p'/T);
     sigma=1e-4;
   svd_trun_index=find(svd_p>sigma^2*1.1);
   [U1,~,~] = svd(p);
 %  ptrunc=U1(:,svd_trun_index)*S1(svd_trun_index,svd_trun_index)*V1(:,svd_trun_index)';
  % Sigma = zf * ptrunc'*pinv(p_correct)/T;
   Sigma=zf*pinv(p);
   Sigmazp = Sigma(:, (m + p1 + p1 + 1):end - 1);
       % [U2,~,~] = svd(p_correct);
        nu=U1(:,max(svd_trun_index)+1:end - 1)';
 %       if size(nu,1)>0
    %    nu = null(p')';
        nuk = nu(:, (m + p1 + p1 + 1):end - 1);
        desired_poles = 0.5 * rand(size(Sigmazp, 1), 1)+0.1;
        


     if ((max(abs(eig(Sigmazp)))) < 1)
         disp('Order of the functional observer');
         disp(size(L,1));
         order= true;
     elseif ~isempty(nuk) && (IsObservable(Sigmazp,nuk))
%         KK = place(Sigmazp', nuk', desired_poles);
        
setlmis([]);
n1=size(zf,1);n2=size(nuk,1);
P = lmivar(1,[n1 1]); % P是一个对称正定矩阵变量
W = lmivar(2,[n1 n2]); % K是一个矩阵变量
Q=Sigmazp;
CLMI1=newlmi;% i=1
lmiterm([CLMI1 1 1 P],-1.1,1); % P的系数为1
lmiterm([CLMI1 1 1 0],Q'*Q); % -Q'Q的系数为-1和1
lmiterm([CLMI1 1 1 -W],nuk',Q,'s'); % -Delta'W'Q - QWDelta
% lmiterm([CLMI1 1 1 W],-Q',nuk); % -Delta'W'Q - QWDelta
lmiterm([CLMI1 1 2 -W],nuk',1); % Delta'W'
lmiterm([CLMI1 2 2 P],-1,1); % P
% 求解LMI
lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
if tmin <0
    P_sol = dec2mat(lmisys,xfeas,P);
    W_sol = dec2mat(lmisys,xfeas,W);
    K_sol = P_sol\W_sol; % 根据K = P^(-1)W计算K
    KK=K_sol';
else
    disp('LMI没有可行解');
end
        
           order= true;
  %      eig(Sigmazp - KK' * nuk);
        Sigma = Sigma - KK' * nu;
   %     Sigmazp = Sigma(:, (m + p1 + p1 + 1):end);
    %    eig(Sigmazp);
        
     disp('Order of the functional observer');
     disp(size(L,1));
     
     else
         disp('The Darouch observer does not exist');
         Sigma=[];
         order=false;
     end

end