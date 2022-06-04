% D^-1 * M D
% ----------
% [      a, (b*y)/x, (c*z)/x]
% [(d*x)/y,       e, (f*z)/y]
% [(g*x)/z, (h*y)/z,       k]

% D = [x y z]
 

%% try 1 - bad result

T = [k2 k1 0; k2 k1 0; 1 0 0];
JJ = inv(I)*simmat(J,T)*I;
u1 = matmis(JJ, 'L1');
%u2 = matmis(JJ, 'L2');
uinf = matmis(JJ, 'Linf');
displayf(JJ, u1, uinf)

%% try 2 - fail

T = [k2 k1 -1; k2 k1 1; 1 0 0];
JJ = inv(I)*simmat(J,T)*I;
u1 = matmis(JJ, 'L1');
%u2 = matmis(JJ, 'L2');
uinf = matmis(JJ, 'Linf');
displayf(JJ, u1, uinf)

%% try 3 - 

T = [k2 k1 -1; 1/(2*k1) 1/(2*k1) -k2/k1; -1/2 1/2 0];
JJ = inv(I)*simmat(J,T)*I;
u1 = matmis(JJ, 'L1');
%u2 = matmis(JJ, 'L2');
uinf = matmis(JJ, 'Linf');
displayf(JJ, u1, uinf)
pretty(JJ)

%% try 4 - chack next

T = [0 0 1; -1 1 0; 1 0 0];
JJ = inv(I)*simmat(J,T)*I;
u1 = matmis(JJ, 'L1');
%u2 = matmis(JJ, 'L2');
uinf = matmis(JJ, 'Linf');
displayf(JJ, u1, uinf)

%% try 4.1 - chack next

T = [0 0 1; -1 1 0; 1 0 0];
DT = diag([d2 k2 1]);
T = T * DT;
JJ = simmat(J,T);
u1 = matmis(JJ, 'L1');
%u2 = matmis(JJ, 'L2');
uinf = matmis(JJ, 'Linf');
displayf(JJ, u1, uinf)

%% try 4.2 - 

T = [0 0 1; -1 1 0; 1 0 0];
DT = diag([-1 1 1]);
T = T * DT;
JJ = simmat(J,T);
u1 = matmis(JJ, 'L1');
%u2 = matmis(JJ, 'L2');
uinf = matmis(JJ, 'Linf');
displayf(JJ, u1, uinf)

%% try 5 - check next

J = [-a 0 -3*k2*asd; 0 -b 3*k2*asd; -1 1 0]; % our matrix
T = [1 0 0; 0 1 0; 0 5 1];
JJ = simmat(J,T);
u1 = matmis(JJ, 'L1');
%u2 = matmis(JJ, 'L2');
uinf = matmis(JJ, 'Linf');
displayf(JJ, u1, uinf)

%% try 5 - check next

T = [0 0 1; 1 1 0; -1 0 0];
JJ = simmat(J,T);
u1 = matmis(JJ, 'L1');
%u2 = matmis(JJ, 'L2');
uinf = matmis(JJ, 'Linf');
displayf(JJ, u1, uinf)




%% T is diagonal matrix

syms t1 t2 t3;
I = diag([1,1,1]);

syms a b c
J=[a 0 -c; 0 b c; -1 1 0];

J_new = inv(T) .* J .* T;

syms d e f g h k
B=[a b c; d e f; g h k];

%%

tryall(J,NaN,NaN)


%%

numzeros(I)
%% 
clc;clear;close all;
syms d1 d2 k1 k2 asd

T = [k2 k2 -1; k2 k1 1; 1 0 0];
J = [-k1/d1 0 -3*k1*asd; 0 -k2/d2 3*k2*asd; -1 1 0];

J1 = inv(T) * J * T;
J2 = T * J * inv(T)
J2*(k1+k2)

%% general
clear all; close all; clc;
%print -depsc

syms a b c d e f g h k x y z
syms d1 d2 k1 k2 asd

M = [a b c; d e f; g h k];
D = diag([x y z]);
I = diag([1 1 1]);
J = [-k1/d1 0 -3*k1*asd; 0 -k2/d2 3*k2*asd; -1 1 0]; % our matrix


% gives the similar matrix of matrix A using transform matrix T
function P = simmat(A, T)
    P = inv(T) * A * T;
end

% gives the matrix measure of matrix A for L1,L2,Linf
function u = matmis(A, L)
    [n m] = size(A);
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
        (A + A.')/2
        u = eig((A+A')/2);
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

% plot formulas formatted as latex
function displayf(varargin)
    figure; hold on; hold off;
    for k = 1:nargin
        str = latex(varargin{k});
        subplot(nargin,1,k); hold all;
        axis off;
        text(0.5, 0.5, ['$$' str '$$'], 'Interpreter','latex', 'FontSize',15, ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end
    grid off; hold off;
end


%from,to in decimal
function x = tryall(J, from, times)
    if (isnan(times))
        times = 19683;
    end
    if (isnan(from))
        from = 0;
    end
    y = from;
    x = 0;
    c = CounterInBase(3);
    c.set(from);
    while (length(c.next(0)) <= 9 && y<=from+times) 
        a = [zeros(1, 9-length(c.next(0))) c.next(0)];
        T = sym(reshape(a,[3,3])) - ones(3,3);
        JJ = simmat(J,T);
        D = diag(JJ);
        if( ~isnan(sum(JJ,'all')) && D(1)~=0 && D(2)~=0 && D(3)~=0 )
%           displayf(T,JJ);\
            if(JJ(1,3)==0 && JJ(3,1)==0)
                x = x+1;
                y
            end
        end
        y = y+1;
%         [x y]
        c.next();
    end
end

% count zero cell in matrix
function n = numzeros(A)
    n = 0;
    for r = 1:size(A,1)
        for c = 1:size(A,1)
            if(A(r,c) == 0)
                n = n+1;
            end
        end
    end
end

