% ---------------------------------------
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
% ---------------------------------------

function [d, a, b] = Triple(K, X)
% K - The number of source symbols
% X - An encoding symbol ID

    global J ;
    global Q ;
    global L_prime ;
    J_K = J(K-4+1);
    A = mod( 53591 + J_K*997 , Q );
    B = mod( 10267 * (J_K+1) , Q );
    Y = mod( B + X*A , Q );
    v = my_rand(Y, 0, 2^20);
    d = Deg(v);
    a = 1 + my_rand(Y, 1, L_prime-1);
    b = my_rand(Y, 2, L_prime);
    
end