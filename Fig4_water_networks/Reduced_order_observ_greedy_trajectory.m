function [Sigma,exist,newL]=Reduced_order_observ_greedy_trajectory(u_traj,X,C,L,~)
n=size(X,1);
up = u_traj(:,1:end-1);
xp = X(:, 1:end-1);
xf = X(:, 2:end);
yp = C * xp;
yf = C * xf;
zp = L * xp;
zf = L * xf;
p = [up; yp; yf; zp];
%f = [p; zf];
E0=eye(n);
yes=data_functional_observability_test(u_traj,C*X,L*X,n);
Sigma=[];
exist=false;
newL=L;
if yes
     disp('This system is functionally observable from data');
     else
     disp('Attension: This system is NOT functionally observable from data');
end
    flag=0;
  % if rank_svd(p) == rank_svd(f)
    if   check_consistent(p,zf)
        [Sigma,yes2] = Darouch_observer_with_trajectory_v3(X,u_traj,C,newL);
        if yes2
            exist=true;
        end
    else
    for i=1:n-size(L,1)
        sumL = sum(L,1);
        r = find(sumL == 0);
        deltaRank = n;
        for j=1: n - size(L,1)
            add_index = r(j);
            Ladd = E0(add_index,:);
            newL=[L;Ladd];
            zp = newL * xp;
            zf = newL * xf;
            p = [up; yp; yf; zp];
            f = [p; zf];
            RP = rank_svd(p);
            RF = rank_svd(f);
            delta = RF - RP;
            if RP==RF
                [Sigma,yes2]=Darouch_observer_with_trajectory(X,u_traj,C,newL);
                if yes2
                    exist=true;
                    flag=1;
                    % newL=[L;Ladd];
                    break
                end
            end
            if delta < deltaRank
                deltaRank = delta;
                L_next = newL;
            end
        end
        if flag==1
          break
        end
           L = L_next;
    end
    if flag==0
        disp('Attension: a reduced order functional observer has not been found');
    end 
    end
end

