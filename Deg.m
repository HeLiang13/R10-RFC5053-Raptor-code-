% ---------------------------------------
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
% ---------------------------------------

function d = Deg(v)
% The degree generator Deg[v] is defined as follows, where v is an
% integer that is at least 0 and less than 2^^20 = 1048576.
% In Table 1, find the index j such that f[j-1] <= v < f[j]
% Then, Deg[v] = d[j]
% +---------+---------+------+
% | Index j | f[j] | d[j] |
% +---------+---------+------+
% | 0 | 0 | -- |
% | 1 | 10241 | 1 |
% | 2 | 491582 | 2 |
% | 3 | 712794 | 3 |
% | 4 | 831695 | 4 |
% | 5 | 948446 | 10 |
% | 6 | 1032189 | 11 |
% | 7 | 1048576 | 40 |
% +---------+---------+------+


    table_1 = [
        0 , 0 , 	0 ;
        1 , 10241 , 1 ;
        2 , 491582 , 2 ;
        3 , 712794 , 3 ;
        4 , 831695 , 4 ;
        5 , 948446 , 10 ;
        6 , 1032189 , 11 ;
        7 , 1048576 , 40 ];
    fj = table_1(:,2);
    dj = table_1(:,3);
    d = [];
    for j = 2:8
        if v < fj(j)
            d = dj(j);
            break;
        end
    end

    if isempty(d)
        error('the input %i is out of range.',v);
    end

end