function Y=Hankel(X,t,N) 

% t is the number of row blocks, N is the number of column blocks
row = size(X(:,1),1);
col = size(X(:,1),2);
Y = zeros(t*row,N*col);
for i = 1:t
    for j = 1:N
        Y((i-1)*row+1:i*row,(j-1)*col+1:j*col) = X(:,i+j-1);
    end
end