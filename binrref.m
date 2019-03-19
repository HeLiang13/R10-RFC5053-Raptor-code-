function [A]=binrref(A,augment_flag)
% binrref��������Ԫ����A���ݻ�,�ı�����ԭʼ��ʵ�����ݻ�
% augment_flag = 1 ��ʾAΪ�������
% Created on Thu Mar 14 18:09:44 2019
% @author: HeL
[m,n] = size(A);
col = n;
if augment_flag
    col = n - 1;
end
i = 1;
j = 1;
while (i <= m) && (j <= col)
    % Find index of first 1 in the remainder of column j.
    k =find((A(i:m,j)),1,'first');
    if ~isempty(k)
        k = k+i-1;
        % Swap i-th and k-th rows.
        A([i k],j:n) = A([k i],j:n);
        % Subtract multiples of the pivot row from all the other rows.
        for k = [1:i-1 i+1:m]
            if(A(k,j))
            A(k,j:n) = bitxor(A(k,j:n),A(i,j:n));
            end
        end
        i = i + 1;
        j = j + 1;
    else
        j=j+1;
    end
  
end

