% ---------------------------------------
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
% ---------------------------------------
function [B] = multiply_bin(A,vec)
% Multiplication on GF(2).
% A is a K*N matrix on GF(2), vec is a N*M matrix.
% B is a K*M matrix.


[K, NA] = size(A);
[Nv, M] = size(vec);

if NA ~= Nv
    error('Inner matrix dimensions must agree.')
end


B = zeros(K, M);

for ii = 1 : K
    for jj = 1 : M
        vec_A = A(ii,:);
        vec_v = vec(:,jj);
        tmp_B = 0;
        index_one = find(vec_A==1);
        for kk = 1 : length(index_one)
            index = index_one(kk);
            tmp_B = bitxor( tmp_B, vec_v(index) );
        end
        B(ii,jj) = tmp_B;
    end    
end


