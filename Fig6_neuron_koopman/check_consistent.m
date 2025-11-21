% =========================================================================
% Filename: check_consistent.m
% Purpose: Check the consistency condition for observer design - verify if zf lies in the row space of p
% Inputs:
%   p  - Data matrix containing input and output information
%   zf - Target functional state matrix to be estimated
% Outputs:
%   exist - Boolean result indicating whether consistency condition is satisfied
% Main functions:
%   Verify the rank condition: rank([p; zf]) == rank(p)
% =========================================================================

function  exist=check_consistent(p,zf)
% check if rank([p;zf])=rank(p); 
resu=zf*(eye(size(zf,2))-pinv(p)*p);
if max(max(abs(resu)))<1e-4
    %  max(max(resu))<1e-3
    exist=true;
else
    exist=false;
end

 