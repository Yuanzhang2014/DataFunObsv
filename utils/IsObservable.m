function flag = IsObservable(A,C)

% This function checks the observability of pair (A,C) when the system
% order is relatively high
% flag = true --> observable
% flag = false --> unobservable

% Author: Ziyuan Luo
% Email: ziyuan.luo@bit.edu.cn
% Last revision: Oct.20, 2025

%------------- BEGIN CODE --------------
    n = size(A,1);
    p = size(C,1);
    Oac = zeros(p*n,n); 
    tmp = C;
    for i = 1:n
        Oac((i-1)*p+1:i*p,:) = tmp;
        tmp = tmp*A;
        tmp = tmp / max(max(tmp));
    end

    if rank(Oac) == n 
        flag = true;
    else
        flag = false;
    end
%------------- END OF CODE --------------

end