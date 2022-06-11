clear all; close all; clc;
%  SYSTEM DEFINITIONS  %

% define system constants
global K D a_21 P_ref ws P_L1 P_L2
K = 8.73e-6; %0.7 works best
D = 0.18 * K;
P_ref = 4e6;
ws = 50 * 2 * pi;
a_21 = (12e3)^2 * abs(-1 / (1j * ws * 0.0275));

% define boundaries
omega_min = 49.8 * 2 * pi; % rad/sec
omega_max = 50.2 * 2 * pi; % rad/sec
delta_min = -pi / 2; % rad
delta_max = pi / 2; % rad

% define resolutions
omega_res = 0.1;
delta_res = 0.1;

%{ 
for 2 gen case
p_ref1 = p_ref;
p_ref2 =p_ref;
P_L1 = p_ref;
P_L2 = p_ref;
K1 =K;
K2 = K;
D1=D;
D2=D;
%}
% F = @(d2, w1, w2, P_1, P_2) [w2-w1 ; K1(3*p_ref1 -3*P_1 - )];


% define Jacobian formula
J = @(d, w) [0, 1; -3 * K * a_21 * cos(d), -K / D]; % d is up

% define numerical grid
[delta_2, omega_2] = meshgrid(delta_min:delta_res:delta_max, omega_min:omega_res:omega_max);

%  ALGORITHM  %

% 1. Solve the power flow for the system (the solution is an equilibrium point of the dynamics).
equilibrium_points = solve_power_flow();
assert(~isempty(equilibrium_points), "The power flow system does not have a solution.")

for i = 1:size(equilibrium_points, 2)
    fprintf("Evaluate equilibrium point (%d,%d) \n", equilibrium_points{i}(1), equilibrium_points{i}(2));
    % 2. Evaluate the Jacobian at the equilibrium point.
    args = num2cell(equilibrium_points{i}); % use the equilibrium point as argument list
    J_eq = J(args{:});

    % 3. Check that the equilibrium point is stable (by checking the eigenvalues of the Jacobian there).
    eigenvalues = eig(J_eq);

    % 4. If it is not stable, stop
    if max(real(eigenvalues)) > 0
        fprintf("The equilibrium point is  not stable.");
        continue;
    end

    % 5. Solve the Lyapunov equation at the equilibrium point using MATLAB’s lyap, using Q = I.
    Q = eye(ndims(J_eq));
    P = lyap(J_eq, Q);

    % 6. Find the square root of the result of lyap using sqrtm. The transformation we will use to define the weighted Euclidean norm is the inverse of the square root.
    T = sqrtm(P);
    A = @(d, w) (T^(-1) * J(d, w) * T);

    % 7. Use meshgrid to define a discrete grid over the state-space, and at each point of the grid calculate the matrix measure induced by the weighted Euclidean norm.
    % TODO make grid_matmis generic for any mesh grid
    %     mu_matrix_L1 = grid_matmis(A, delta_2, omega_2, 'L1');
    mu_matrix_L2 = grid_matmis(A, delta_2, omega_2, 'L2');
    %     mu_matrix_Linf = grid_matmis(A, delta_2, omega_2, 'Linf');

    % 8. Use contour (or any other function of your choice) to draw the boundary of the convergence, that is, the line along which the matrix measure is zero.
    % TODO calculate distances function and convert it back to the original norm
    % TODO plot for all mu
    dist = sqrt(delta_2.^2 + omega_2.^2);
    plot_contraction(delta_2, omega_2, mu_matrix_L2, dist, equilibrium_points{i});
    figure;
    hold on;
    plot(delta_2, mu_matrix_L2(1, :));
    hold off;
end

% for base example, don't need to solve it just return the known points
function equilibrium_points = solve_power_flow()
    global P_ref a_21
    equilibrium_points = cell([1, 2]);
    equilibrium_points{1} = [(asin(P_ref / a_21)), 0]; % (delta,omega) point
    equilibrium_points{2} = [pi - (asin(P_ref / a_21)), 0]; % (delta,omega) point
end

%{
function equilibrium_points = solve_power_flow(delta_2, omega_2)
    global P_ref a_21 P_L1 P_L2
    P_1 = -a_21 * sin(delta_2) +P_L1;
    P_2 = a_21 * sin(delta_2) +P_L2;

end
%}