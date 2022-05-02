

delta = 0:180;
omega = 48:50;

mu_matrix = mu_lyaponov_one_gen (delta,K,D,a_12);

%show mu_matrix in 3D
%show areas when mu_matrix<0

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
    




