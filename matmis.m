
% gives the matrix measure of matrix A for L1,L2,Linf
function mu = matmis(A, L)
    [row_num, col_num] = size(A);
    
    % L1 -> take the 'max' column  
    if strcmp(L, 'L1')
        col_max = NaN;
        for col = 1:col_num
            col_sum = sum(abs(A(1:col-1,col))) + A(col,col) + sum(abs(A(col+1:col_num,col)));
            if isnan(col_max) || col_sum > col_max
                col_max = col_sum;
            end
        end
        mu = col_max;
        return;
    end
    
    % Linf-> take the 'max' row 
    if strcmp(L,'Linf')
        row_max = NaN;
        for row = 1:row_num
            row_sum = sum(abs(A(row,1:row-1))) + A(row,row) + sum(abs(A(row,row+1:row_num)));
            if isnan(row_max) || row_sum > row_max
                row_max = row_sum;
            end
        end
        mu = row_max;
        return;
    end
    
    % L2 -> take the 'max' eigenvalue
    if strcmp(L, 'L2')
        mu = max(eig((A + A.')/2));
        return;
    end
    
    assert(false, 'L aregument is wrong');  
end

