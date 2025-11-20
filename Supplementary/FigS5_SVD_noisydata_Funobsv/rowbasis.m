function T=rowbasis(M) % find the row basis of a matrix using SVD
[~,~,V] = svd(M);
 r=rank_svd(M);
 T=V(:,1:r)';
end