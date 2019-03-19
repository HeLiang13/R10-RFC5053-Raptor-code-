% ---------------------------------------
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
% ---------------------------------------

function [C_solve] = solve_by_decoder(A_mat, D)
% Equation is D = A_mat * C, on GF(2)
% To solve C = A_mat^(-1) * D
% Decoder of Raptor code in RFC5053 is applied.



method_index = 1;
decoding_methods = { 'standard' , 'Kim' , 'random' };
decoding_method = decoding_methods{method_index};

[M, L] = size(A_mat);

% calculate row_score for Kim's algorithm
if strcmp(decoding_method,'Kim')
    row_score = zeros(M,1);
    for jj = 1 : M
        logic_index = A_mat(jj,:) == 1;
        tmp_mat = A_mat(:,logic_index);
        row_score(jj) = sum(sum(tmp_mat));
    end    
end

%%
%%%%%%%%%%%%---------First Phase ---------%%%%%%%%%%%  
i = 0;
u = 0;
row_degree_of_V_ori = sum(A_mat,2);     % row_degree_of_V(_ori) is a vector of length M-i
row_degree_of_V = row_degree_of_V_ori;  % for original degree, U is considered. 
col_order = 1 : L;


while (1)
    % calculate r and index_min_deg
    r = min(row_degree_of_V);
    index_min_deg = find(row_degree_of_V == r);  % notice that index_min_deg is indexed form 1 to M-i
    if r == 0                        % r should be positive unless V is of all zeros.
        tmp_vec = row_degree_of_V;
        tmp_ind = 1 : length(row_degree_of_V);
        tmp_vec(index_min_deg) = [];
        tmp_ind(index_min_deg) = [];
        if ~isempty(tmp_vec)
            r = min(tmp_vec);
            index_min_deg = tmp_ind(tmp_vec == r);            
        else
            break;
        end
    elseif isempty(r)
        break;
    end
    
    % choose a base row
    switch decoding_method
        case 'standard'
            if r ~= 2
                tmp_degree = row_degree_of_V_ori( index_min_deg );
                [~,tmp_index] = min(tmp_degree);
                chosen_row_index = index_min_deg(tmp_index) ;
            elseif r==2
                Mat_d2_inV = A_mat(i+index_min_deg , i+1:L-u);
                comp_row_index = find_component(Mat_d2_inV);
                index_min_deg = index_min_deg(comp_row_index);
                chosen_row_index = index_min_deg(1);
            end
        case 'Kim'
            if r == 1
                tmp_rand = randi(length(index_min_deg));
                chosen_row_index = index_min_deg(tmp_rand) ;
            else
                tmp_score = row_score(index_min_deg + i);
                [~,tmp_index] = max(tmp_score);
                chosen_row_index = index_min_deg(tmp_index);
            end
        case 'random'
            tmp_rand = randi(length(index_min_deg));
            chosen_row_index = index_min_deg(tmp_rand) ;
    end
    
    % swap rows
    A_mat([i+1,i+chosen_row_index],i+1:L) = A_mat([i+chosen_row_index,i+1],i+1:L); 
    D([i+1,i+chosen_row_index]) = D([i+chosen_row_index,i+1]);
    row_degree_of_V_ori([1,chosen_row_index]) = row_degree_of_V_ori([chosen_row_index,1]);
    row_degree_of_V([1,chosen_row_index]) = row_degree_of_V([chosen_row_index,1]);
    first_row = A_mat(i+1,i+1:L-u);
    if strcmp(decoding_method,'Kim')
        row_score([i+1,i+chosen_row_index]) = row_score([i+chosen_row_index,i+1]);
    end
    
    % swap first columns
    ex_col_index = find(first_row);
    ex_col_index_first = ex_col_index(1);
    A_mat( i+1:M , [i+1,i+ex_col_index_first] ) = A_mat( i+1:M , [i+ex_col_index_first,i+1] );
    col_order([i+1,i+ex_col_index_first]) = col_order([i+ex_col_index_first,i+1]);    
    
    % swap remaining columns    
    rem = r - 1;
    if rem > 0
        for kk = 1:rem
            tmp_ex_col_index = ex_col_index(end-kk+1);
            A_mat( i+1:M ,[L-u-kk+1,i+tmp_ex_col_index] ) = A_mat( i+1:M ,[i+tmp_ex_col_index,L-u-kk+1] );
            col_order([L-u-kk+1,i+tmp_ex_col_index]) = col_order([i+tmp_ex_col_index,L-u-kk+1]);
        end        
    end
    
    % update relevant varibles
    i = i + 1;
    u = u + rem;
    row_degree_of_V_ori(1) = [];
    row_degree_of_V(1) = [];
    if rem > 0 
        last_columns = A_mat(i+1:M , L-u+1:L-u+rem);
        minus_degree = sum(last_columns,2);
        row_degree_of_V = row_degree_of_V - minus_degree;
    end
    
    % xor operation
    base_row = A_mat(i, i : L);
    first_column = A_mat(i + 1 : M , i);
    xor_row_index = find(first_column);
    xor_row_index = xor_row_index + i;
    for jj = 1 : length(xor_row_index)        
        kk = xor_row_index(jj);
        A_mat(kk,i : L) = bitxor( A_mat(kk,i : L), base_row );
        D(kk) = bitxor( D(kk), D(i) );
        tmp_degree_ori = sum(A_mat(kk,i : L));
        tmp_degree = tmp_degree_ori - sum( A_mat(kk,L-u+1:L) );
        row_degree_of_V_ori(kk-i) = tmp_degree_ori;
        row_degree_of_V(kk-i) = tmp_degree;      
    end  
    
end

if i + u ~= L
    error('First Phase fails');
end

%%
%%%%%%%%%%%%---------Second Phase ---------%%%%%%%%%%%  

U_lower = A_mat(i+1:M,i+1:L);
D_lower = D(i+1:M);
result = binrref( [U_lower , D_lower], 1 );
U_lower = result(1:u,1:u);
D_lower = result(1:u,u+1);

if ~isequal(U_lower , eye(u))
    error('Second Phase fails');
end

A_mat(i+1:L,i+1:L) = U_lower;
A_mat(L+1:M,:) = [];
D(i+1:L) = D_lower;
D(L+1:M) = [];

%%
%%%%%%%%%%%%---------Third/Fourth Phase ---------%%%%%%%%%%%  

U_upper = A_mat(1:i,i+1:L);
D_upper_xor = multiply_bin(U_upper,D_lower);
D(1:i) = bitxor( D(1:i) ,D_upper_xor);
A_mat(1:i,i+1:L) = 0;

%%
%%%%%%%%%%%%--------- Solving Intermediate Symbols and Source Symbols ---------%%%%%%%%%%%  

% Column swaps have been employed to transfer A_mat into Identity Matrix,
% which is recorded in variable "col_order". Then we have 
% A_mat * C(col_order) = D, that is C(col_order) = D.

C_solve = zeros(L,1);
for jj = 1 : L
    C_solve(jj) = D(col_order == jj);
end




end