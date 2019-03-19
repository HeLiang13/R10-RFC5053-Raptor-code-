% ---------------------------------------
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
% ---------------------------------------

function rand_num = my_rand(X,i,m)
% X is a non-negative integer, i is a non-negative integer, and m
% is a positive integer, and the value produced is an integer between 0
% and m-1.


% rand_num = (   V0( mod((X + i),256) ) xor V1( mod((floor(X/256)+ i),256) )   ) mod m
    global V0 ;
    global V1 ;
    input1 = V0( mod(X + i,256) + 1 ); 
    input2 = V1( mod(floor(X/256)+ i,256) + 1);
    temp = bitxor(input1,input2);
    output = mod( temp,m );
    rand_num = double(output);
    
end