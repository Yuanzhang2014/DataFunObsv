function r=rank_svd(M)
% for i=1:size(M,1) % normalize
%     M(i,:)=M(i,:)/(max(M(i,:)));
% end
epsion=1e-6;
S=svd(M);
indx=find(S>epsion);
r=length(indx);
end

