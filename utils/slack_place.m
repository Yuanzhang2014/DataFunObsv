function [K,nu] = slack_place(Sigma,p,T)
    
    thresh = 1e-20;
    Sigmazp = Sigma(:, end-size(Sigma,1)+1:end);
    svd_p= svd(p*p'/T);
    
    e = eig(Sigmazp);
    idx = abs(e) > 1;
    desired_poles = e;
    desired_poles(idx) = .8*rand(sum(idx),1) + 0.1;
    
    % svd_trun_index=find(svd_p>thresh); % BA
    % svd_trun_index=find(svd_p>tol); % auto
    % svd_trun_index=find(svd_p>5e-9); % ER
    % svd_trun_index=find(svd_p>5e-8); % Example
    [U1,~,~] = svd(p);
    % nu=U1(:,max(svd_trun_index)+1:end)';
    % nuk = nu(:, end-size(Sigmazp,1)+1:end);
    warning('off','all');
    for i = max(find(svd_p>thresh)) :-2:1
        nu = U1(:,i+1:end)';
        nuk = nu(:, end-size(Sigmazp,1)+1:end);
        try 
            prec=[];
            [K, prec] = place(Sigmazp', nuk', desired_poles);
        catch ME
        end
        if prec > 0
            break
        end
    end
    warning('on','all');
    % error = max(max(abs(K'*nu*p)));
    % KK = place(Sigmazp', nuk', desired_poles);
    % rank_svd(ctrb(Sigmazp', nuk'))
    % test = Sigma - KK'*nu;
    % max(abs(eig(test(:, end-size(Sigma,1)+1:end))))
    % max(max(abs(nu*p)))
    % max(max(abs(KK'*nu*p)))
end