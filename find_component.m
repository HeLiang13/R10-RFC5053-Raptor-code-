% ---------------------------------------
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
% ---------------------------------------

function    [row_index] = find_component(Mat)
% Mat is a matrix on GF(2) where each row is of degree 2.
% This function output the logical row index of Mat,
% where Mat(row_index,:) is the maximum connected submatrix of Mat. 
% But actually, the degree of Mat can be of any integer.


Component = cell(1,2);  
number_of_component = 1;
[m, n] = size(Mat);
if find( sum(Mat,2) ~= 2 )
    error('Some rows of the input dont have degree of 2.')
end

Mat = Mat == 1;              % transfer Mat into a logical one for later logical indexing.
row_index_done = false(m,1); % indicate the row_index which have been covered.
row_index_todo = false(m,1); % indicate the row_index which will be covered in the next iteration.
row_index_todo(1) = 1;
col_index_done = false(1,n);
row_index = false(m,1);

while(1)
        
    col_index_todo = Mat(row_index_todo,:) ;
    col_index_todo = sum(col_index_todo,1) == 1;
    col_index_todo = xor(col_index_todo,col_index_done);  
    
    [tmp_row_index,~] = find( Mat(:,col_index_todo) );    
    col_index_done(col_index_todo) = 1;
    row_index(tmp_row_index) = 1;
    
    row_index_done(row_index_todo) = 1;
    row_index_todo = xor(row_index , row_index_done);    
    
    if isempty(find(row_index_todo, 1))
        if length(find(col_index_done)) >= n/2  % satisfying this, the output must be the maximum component
            Component{number_of_component,1} = row_index;
            Component{number_of_component,2} = length(find(col_index_done));
            break;
        else
            % store the component for comparison at the end
            Component{number_of_component,1} = row_index;
            Component{number_of_component,2} = length(find(col_index_done));     
            % reset variables
            row_index_done = false(m,1); 
            row_index_todo = true(m,1);
            col_index_done = false(1,n);
            row_index = false(m,1);
            % skip the row_index that alrealy exist in former component
            for ii = 1 : number_of_component
                tmp_row_index = Component{ii,1};
                row_index_todo(tmp_row_index) = 0;
            end           
            
            tmp_ind = find(row_index_todo,1);
            if isempty(tmp_ind)
                break;
            end
            
            row_index_todo = false(m,1); 
            row_index_todo(tmp_ind) = 1;            
            number_of_component = number_of_component + 1;
        end
    end
    
end

if number_of_component >=2
    size_component = zeros(1,number_of_component);
    for ii = 1 : number_of_component
        size_component(ii) = Component{ii,2};
    end
    [~,tmp_ind] = max(size_component);
    row_index = Component{tmp_ind,1};    
end



% end