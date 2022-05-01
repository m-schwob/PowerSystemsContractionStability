classdef CalcTool
    methods (Static)
        % gives the similar matrix of matrix A using transform matrix T
        function P = simmat(A, T)
            P = inv(T) * A * T;
        end

        % gives the matrix measure of matrix A for L1,L2,Linf
        function u = matmis(A, L)
            [n, m] = size(A);
            assert(n == m, "matrix isn't squre");
            syms none

            if strcmp(L, 'L1')
                rows_sum = ones(n,1) * none;
                for c = 1:n
                    rows_sum(c) =  A(c,c);
                    for r = 1:n
                        if r ~= c   
                            rows_sum(c) = rows_sum(c) + abs(A(r,c));
                        end
                    end
                end
                u = rows_sum;
                return;
            end

            if strcmp(L, 'L2')
                u = (A + A.')/2;
                %u = eig((A+A')/2);
                return;
            end

            if strcmp(L,'Linf')
                rows_sum = ones(n,1) * none;
                for r = 1:n
                    rows_sum(r) =  A(r,r);
                    for c = 1:n
                        if r ~= c   
                            rows_sum(r) = rows_sum(r) + abs(A(r,c));
                        end
                    end
                end
                u = rows_sum;
                return;
            end

            assert(false, 'L aregument is wrong');  
        end
    end
end
