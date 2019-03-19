% Decoder of Raptor code in RFC5053( 3GPP Multimedia Broadcast/Multicast Services standard)
% ---------------------------------------
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
% ---------------------------------------

% function   [suc_flag, time_pass, u  ] = decode (K, G_LT_all, Symbol_all, method_index)

global G_LDPC;
global G_Half;
global S;
global H;
global L;
global Num_repair;


% VideoObject = VideoWriter('DecodingPhase.avi');  % The decoding process is recorded and made into Video
% VideoObject.FrameRate = 100 ;
% open(VideoObject);


% choose decoding method for First Phase.
% when r=2, standard algorithm choose row in the maximum size component.
% Kim's method is in paper of DOI 10.1109/LCOMM.2008.080599;
% Kim's choose row when r>=2 with the highest score.
% random method choose row at random.
% Zhang's method is in paper of DOI 10.1109/ICMTMA.2010.531.
% A_mat should be a M*L matrix.


method_index = 1;
decoding_methods = { 'standard' , 'Kim' , 'random', 'Zhang' };
decoding_method = decoding_methods{method_index};

% Choose receive symbols.
% Let N >= K be the number of received encoding symbols and then M = S+H+N.
% Let C be the column vector of the L intermediate symbols, 
% and let D be the column vector of M symbols.

overhead = 1.05;
N = ceil(overhead * K);
Received_Symbol_Index = randperm(K+Num_repair, N);
Received_Symbol = Symbol_all(Received_Symbol_Index);  % Received_Symbol is a row vector
G_LT_receive = G_LT_all(Received_Symbol_Index,:);
I_S = eye(S);
I_H = eye(H);
A_mat = [
        G_LDPC, I_S, zeros(S,H);
        G_Half,  I_H;
        G_LT_receive
        ];
N = length(Received_Symbol);
M = size(A_mat,1);  
D = [zeros(1,H + S), Received_Symbol]';

% calculate row_score for Kim's algorithm
if strcmp(decoding_method,'Kim')
    row_score = zeros(M,1);
    for jj = 1 : M
        logic_index = A_mat(jj,:) == 1;
        tmp_mat = A_mat(:,logic_index);
        row_score(jj) = sum(sum(tmp_mat));
    end    
end

time1 = clock;
%%
%%%%%%%%%%%%---------First Phase ---------%%%%%%%%%%%  
i = 0;
u = 0;
row_degree_of_V_ori = sum(A_mat,2);     % row_degree_of_V(_ori) is a vector of length M-i
row_degree_of_V = row_degree_of_V_ori;  % for original degree, U is considered. 
col_order = 1 : L;


%  preprocessing for Zhang's algorithm
if strcmp(decoding_method,'Zhang')    
    A_mat(:,[1:S+H,K+1:L]) = A_mat(:,[K+1:L,1:S+H]);
    col_order([1:S+H,K+1:L]) = col_order([K+1:L,1:S+H]);     
    i = S + H;
    
    % zero out submatrix under I_S, I_H
    mul_mat1 = A_mat(S+1:M, 1:S);
    mul_mat2 = A_mat(S+H+1:M, S+1:S+H);
    A_mat(S+1:M,S+H+1:L) = bitxor( A_mat(S+1:M,S+H+1:L) , multiply_bin(mul_mat1, A_mat(1:S,S+H+1:L)) );
    D(S+1:M) = bitxor( D(S+1:M) , multiply_bin(mul_mat1, D(1:S)) );
    A_mat(S+1:M, 1:S) = 0;
    A_mat(S+H+1:M,S+H+1:L) = bitxor( A_mat(S+H+1:M,S+H+1:L) , multiply_bin(mul_mat2, A_mat(S+1:S+H,S+H+1:L)) );
    D(S+H+1:M) = bitxor( D(S+H+1:M) , multiply_bin(mul_mat2,D(S+1:S+H)) );
    A_mat(S+H+1:M, S+1:S+H) = 0;
    
    row_degree_of_V_ori = sum(A_mat(i+1:M,:),2);
    row_degree_of_V = row_degree_of_V_ori ; 
end

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
    % When r=2, RFC choose row with exactly 2 ones in V that is
    % part of a maximum size component in the graph defined by V.
    % This is where the  computational bottleneck lies.
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
        case {'random','Zhang'}
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
    A_mat( : , [i+1,i+ex_col_index_first] ) = A_mat( : , [i+ex_col_index_first,i+1] );
    col_order([i+1,i+ex_col_index_first]) = col_order([i+ex_col_index_first,i+1]);    
    
    % swap remaining columns    
    rem = r - 1;
    if rem > 0
        for kk = 1:rem
            tmp_ex_col_index = ex_col_index(end-kk+1);
            A_mat( : ,[L-u-kk+1,i+tmp_ex_col_index] ) = A_mat( : ,[i+tmp_ex_col_index,L-u-kk+1] );
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
    
    % get the current frame and write it into the video. %%%%%%% aditional code 
%     imshow(1-A_mat);
%     frame = getframe(gcf);
%     writeVideo(VideoObject,frame);
    
end

% if i + u == L
%     disp('First Phase succeeds');
% else
%     error('First Phase fails');
% end

% figure(1);
% imshow(1-A_mat);

%%
%%%%%%%%%%%%---------Second Phase ---------%%%%%%%%%%%  

% The submatrix U is further partitioned into the first i rows,
% U_upper, and the remaining M - i rows, U_lower.
% Gaussian elimination is performed on U_lower.
% Convert it into a the identity matrix on first u rows. 

U_lower = A_mat(i+1:M,i+1:L);
D_lower = D(i+1:M);
result = binrref( [U_lower , D_lower], 1 );
U_lower = result(1:u,1:u);
D_lower = result(1:u,u+1);

% if isequal(U_lower , eye(u))
%     disp('Second Phase succeeds');
% else
%     error('Second Phase fails');
% end

A_mat(i+1:L,i+1:L) = U_lower;
A_mat(L+1:M,:) = [];
D(i+1:L) = D_lower;
D(L+1:M) = [];


% figure(2);
% imshow(1-A_mat);

% get the current frame and write it into the video.%%%%%%% aditional code 
% [Height,Width,~] = size(frame.cdata);
% imshow(1-A_mat);
% frame = getframe(gcf);
% frame.cdata = imresize(frame.cdata, [Height Width]);  % keep the frame size the same as the former, if it differs, error may occur
% writeVideo(VideoObject,frame);



%%
%%%%%%%%%%%%---------Third/Fourth Phase ---------%%%%%%%%%%%  

% The method introduced in RFC5053 for Third/Fourth Phase is somehow
% cumbersome. It seems that, zeroing out U_upper directly by U_lower is
% equivalent to the introduced submatrix-like methods, and it's convenient.
% Actually,after the zeroing out, A = L by L identity matrix,
% and D(1:i) = D(1:i) + U_upper * D_lower which is performed on GF(2).

if   strcmp(decoding_method,'Zhang')  
    base_left = A_mat(1:S+H,S+H+1:L-u);
    base_right = A_mat(1:i,L-u+1:L);
    D(1:i) = bitxor( D(1:i) , multiply_bin(base_right,D_lower) );
    D(1:S+H) = bitxor( D(1:S+H) , multiply_bin(base_left,D(S+H+1:L-u)) );    
    A_mat(1:S+H,S+H+1:L-u) = 0;
    A_mat(1:i,L-u+1:L) = 0;
else
    U_upper = A_mat(1:i,i+1:L);
    D_upper_xor = multiply_bin(U_upper,D_lower);
    D(1:i) = bitxor( D(1:i) ,D_upper_xor);
    A_mat(1:i,i+1:L) = 0;    
end

% get the last frame and write it into the video.%%%%%%% aditional code 
% imshow(1-A_mat);
% frame = getframe(gcf);
% frame.cdata = imresize(frame.cdata, [Height Width]);
% writeVideo(VideoObject,frame);
% close(VideoObject);


%%
%%%%%%%%%%%%--------- Solving Intermediate Symbols and Source Symbols ---------%%%%%%%%%%%  

% Column swaps have been employed to transfer A_mat into Identity Matrix,
% which is recorded in variable "col_order". Then we have 
% A_mat * C(col_order) = D, that is C(col_order) = D.

C_solve = zeros(L,1);
for jj = 1 : L
    C_solve(jj) = D(col_order == jj);
end

G_LT_source = G_LT_all(1 : K,:);
Source_Symbol = Symbol_all(1 : K);
Source_Symbol_solve = multiply_bin( G_LT_source , C_solve );

time2 = clock;
time_pass = etime(time2,time1);

if isequal(Source_Symbol_solve',Source_Symbol)
    suc_flag = 1;
else
    suc_flag = 0;
end





