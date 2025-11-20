function newL=compact_firstrowL(oldL,L) % finding newL which is the row basis of oldL, and the first rows of newL is L
% COMPACT_FIRSTROWL - Compute a row basis with specified first rows
% Constructs a new matrix that forms a row basis of oldL with the first rows
% exactly equal to L, ensuring linear independence and proper structure.

% Author: Yuan Zhang
% email: zhangyuan14@bit.edu.cn 
% Last revision: Oct-20-2025
oldL=rowbasis(oldL);
T1=L*pinv(oldL);
T1_orth=null_threshold(T1)';
newL=[L;T1_orth*oldL];
