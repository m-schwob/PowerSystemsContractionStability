close all;
clc;

K = 8.73e-7;
D = 0.18*K;
P_ref = 4e6;
ws = 50*2*pi;
a_12 = (12e3)^2 * abs(-1/(1j*ws*0.0275));
w_res = 0.1;
d_res = 0.1;
delta = linspace(-pi/2,pi/2,1000);
omega = linspace(49.8,50.2,5);
N=length(delta);
mu_matrix = zeros(N,3);
%
for row=1:N
% mu_matrix(row,1) = mu_lyaponov_one_gen (delta(row),K,D,a_12,'L1');
mu_matrix(row,2) = mu_lyaponov_one_gen (delta(row),K,D,a_12,'L2');
% mu_matrix(row,3) = mu_lyaponov_one_gen (delta(row),K,D,a_12,'Linf');
end
mu_matrix = real(mu_matrix);

%show areas when mu_matrix<0
figure (1)
hold on
plot (delta,mu_matrix(:,1))
plot (delta,mu_matrix(:,2))
plot (delta,mu_matrix(:,3))
legend ('L1','L2','Linf')

mu_matrix_cozen = mu_matrix(:,2)';
mu_matrix_cozen = repmat(mu_matrix_cozen,length(omega),1);
figure(2)
hold on
contourf (delta,omega,-mu_matrix_cozen,[0,0])
xlabel('delta')
ylabel('omega')

%{
[x,y] = meshgrid(-5:0.1:5);
f1 = x.^2 + y.^2-10;
f2 = (x+3).^2 + y.^2-10;

hold on;
contourf(x,y,-f1,[0,0]);
contourf(x,y,-f2,[0,0]);
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off;
%}

function mu_matrix = mu_lyaponov_one_gen (delta,K,D,a_12,L)
%     J = [-K/D -cos(delta)*3*a_12*K;
%         1 0];
    J = [0, 1; -cos(2.423481e-01)*3*a_12*K,-K/D];
    Q = eye(2);
    P = lyap (J,Q);
    je =eig(J);
    e = eig(P);
    sqrt_P = sqrtm(P);
    J = [0, 1; -cos(delta)*3*a_12*K,-K/D];
    A = sqrt_P^(-1)*J*sqrt_P;
    mu_matrix = matmis (A,L);
end

% gives the matrix measure of matrix A for L1,L2,Linf
function u = matmis(A, L)
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
        u = col_max;
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
        u = row_max;
        return;
    end
    
    % L2 -> take the 'max' eigenvalue
    if strcmp(L, 'L2')
        u = max(eig((A + A.')/2));
        return;
    end
    
    assert(false, 'L aregument is wrong');  
end





 
