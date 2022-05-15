

delta = 0:180;
omega = 48:50;

mu_matrix = mu_lyaponov_one_gen (delta,K,D,a_12);

%show mu_matrix in 3D
%show areas when mu_matrix<0


%% drawing sample

[x,y] = meshgrid(-5:0.1:5);

f1 = x.^2 + y.^2-10;
f2 = (x+3).^2 + y.^2-10;

%surf(x,y,-f1);

hold on;
C1 = contourf(x,y,-f1,[0,0]);
C2 = contourf(x,y,-f2,[0,0]);
scatter(0,0);
contour(x,y,sqrt(x.^2 + y.^2),[1,1]);
axis square;
xlabel('X')
ylabel('Y')
zlabel('Z')
hold off;
    

%%

function mu_matrix = mu_lyaponov_one_gen (delta,K,D,a_12)
    J = [-K/D -cos(delta)*3*a_12*K;
        1 0];
    Q = eye(2);
    P = lyap (J,Q);
    sqrt_P = sqrtm(P);
    A = sqrt_P*J;
    mu_matrix_L1 = matmis (A,'L1');
    mu_matrix_L2 = matmis (A,'L2');
    mu_matrix_Linf = matmis (A,'Linf');
    mu_matrix = min(mu_matrix_L1,mu_matrix_L2,mu_matrix_Linf);
end



