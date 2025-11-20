function r=rank_svd(M,tol)
% for i=1:size(M,1) % normalize
%     M(i,:)=M(i,:)/(max(M(i,:)));
% end
if nargin == 1
    epsilon = 1e-6;
else
    epsilon = tol;
end
S=svd(M);
indx=find(S>epsilon);
r=length(indx);
end

