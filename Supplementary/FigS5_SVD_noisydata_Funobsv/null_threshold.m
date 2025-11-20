function H=null_threshold(M,tol) % Finding the null space using SVD with threshold 1e-4 (can be adjusted) 
if nargin == 1
    epsilon=1e-4;
else
    epsilon = tol;
end
[~,S,V] = svd(M);
r=size(M,2);
svd_index=find(S>epsilon);

%% rank diff
rr = rank(M);
if rr - rank_svd(M)  ~= 0
    ;
end

%%

H=V(:,size(svd_index,1)+1:r);
%H=H';