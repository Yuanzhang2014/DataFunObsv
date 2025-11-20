function  exist=check_consistent(p,zf)
% Check if rank([p;zf])=rank(p); the default threshold is set to be 1e-7
% but can be adjusted. 
resu=zf*(eye(size(zf,2))-pinv(p)*p);
if max(max(abs(resu)))<1e-4
    %  max(max(resu))<1e-3
    exist=true;
else
    exist=false;
end