function yes=data_functional_observability_test(u_traj,y_traj,z_traj,n) % check the functional observability using data
% n is the system dimension, an upper bound of the observability index 
yes=false;
T=size(y_traj,2);
N=T+1-n; % N is the number of columns
HankU=Hankel(u_traj,n,N);
HankY=Hankel(y_traj,n,N);
HankZ=Hankel(z_traj,n,N);
if rank([HankU;HankY;HankZ])==rank([HankU;HankY])
    yes=true;
end
end

% function Y=Hankel(X,t,N) % t is the number of row blocks, N is the number of column blocks
% Y=[];
% for i=1:t
%     Yi=[];
%     for j=1:N
%         Yi=[Yi,X(:,i+j-1)];
%     end
%     Y=[Y;Yi];
% end
% end
function Y=Hankel(X,t,N) % t is the number of row blocks, N is the number of column blocks
p = size(X,1);
Y = zeros(t * p,N);
for i = 1 : N
    for j = 1 : t
        Y((j - 1) * p + 1 : j * p, i) = X(:, i + j - 1);
    end
end
end


