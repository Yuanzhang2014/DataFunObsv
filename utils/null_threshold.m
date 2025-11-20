function H=null_threshold(M,tol) 
% finding the null space using SVD with threshold 1e-6; 

if nargin == 1
    epsilon=1e-6;
else
    epsilon = tol;
end
[~,S,V] = svd(M);
r=size(M,2);
svd_index=find(S>epsilon);
H=V(:,size(svd_index,1)+1:r);