function [M1,zf,newF]=F_augmention_v3(M1,zf,xp,xf) %% Augmenting F by subspace intersection in each iteration, M1=[Up',Yp',Yf',Zp']'
%%%% finding newF increasing the rank_svd of [M1;newF*xp] but not [M1;zf;newF*xp;newF*xf]
newF=[];
M1=rowbasis(M1);
%if rank_svd(M1)<rank_svd([M1;zf])
if check_consistent(M1,zf)==false   
    n=size(xp,1);
    T1=null_threshold([xp;rowbasis([M1;zf])]');
    T1=T1';
    Phip=T1(:,1:n);
    T2=null_threshold([xf;rowbasis([M1;zf])]');
    T2=T2';
    Phif=T2(:,1:n);
 %   Phir=subspace_union(Phip',Phif');
    
     Phir= SVD_subspace_intersection(Phip',Phif');
    Phir=Phir';
   % T3=null_threshold([Phir*xp;rowbasis(M1)]');
    T3=null_threshold([xp;rowbasis(M1)]');
    T3=T3';
  %  n1=size(Phir,1);
    Phil=T3(:,1:n);
%    Fadd=subspace_union(Phip',Ph)
  %  Fadd=subspace_union(null_threshold(Phil),Phir');
 Fadd= SVD_subspace_intersection(null_threshold(Phil),Phir');
    if ~isempty(Fadd)
        ;
    else 
        ;
    end
    %%update data
    newF=[newF;Fadd'];
    M1=[M1;newF*xp];
    M1=rowbasis(M1);
    zf=rowbasis([zf;Fadd'*xf]);
end
end