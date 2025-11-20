function  exist=check_consistent(p,zf)
% check if rank([p;zf])=rank(p); the threshold is set to be 5e-7
% exist = 1 --> true
% exist = 0 --> false

% Author: Ziyuan Luo
% email: ziyuan.luo@bit.edu.cn 
% Last revision: Oct-20-2025

resu=zf*(eye(size(zf,2))-pinv(p)*p);
if sqrt(trace(resu*resu') / trace(zf*zf')) < 5e-7
    exist=true;
else
    exist=false;
end
