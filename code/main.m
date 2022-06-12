clear all; close all; clc;
%  SYSTEM DEFINITIONS  %

% define system constants
global K D a_21 P_ref ws P_L
K = 8.73e-8; %0.7 works best
D = 0.18 * K;
P_ref = 4e6;
P_L = -5e6;
ws = 50 * 2 * pi;
a_21 = (12e3)^2 * abs(-1 / (1j * ws * 0.0275));


% define boundaries
w_min = 49 * 2 * pi; % rad/sec
w_max = 51 * 2 * pi; % rad/sec
delta_min = -pi / 2; % rad
delta_max = pi / 2; % rad

%simulation parameters
t_final = 10;
t_max_step = 0.01;

% define resolutions
w_res = 0.01;
delta_res = 0.01;

% define Jacobian formula
J = @(d) [0, 1; -3 * K * a_21 * cos(d), -K / D]; % d is up

% define numerical grid
[delta_2, omega_2] = meshgrid(delta_min:delta_res:delta_max, w_min:w_res:w_max);

%  ALGORITHM  %

% 1. Solve the power flow for the system (the solution is an equilibrium point of the dynamics).
equilibrium_points = solve_power_flow();
assert(~isempty(equilibrium_points), "The power flow system does not have a solution.")

for i = 1:size(equilibrium_points, 2)
    fprintf("Evaluate equilibrium point (%d,%d) \n", equilibrium_points{i}(1), equilibrium_points{i}(2));
    % 2. Evaluate the Jacobian at the equilibrium point.
    J_eq = J(equilibrium_points{i}(1));

    % 3. Check that the equilibrium point is stable (by checking the eigenvalues of the Jacobian there).
    eigenvalues = eig(J_eq);

    % 4. If it is not stable, stop
    if max(real(eigenvalues)) > 0
        fprintf("The equilibrium point is  not stable.\n");
        continue;
    end

    % 5. Solve the Lyapunov equation at the equilibrium point using MATLAB’s lyap, using Q = I.
    d_eq = equilibrium_points{i}(1);
    w_eq = equilibrium_points{i}(2);
    Q = eye(ndims(J_eq));
    P = lyap(J_eq, Q);

    % 6. Find the square root of the result of lyap using sqrtm. The transformation we will use to define the weighted Euclidean norm is the inverse of the square root.
    T = sqrtm(P);
    P_final = T^-1;
    A = @(d, w) (T^(-1) * J(d) * T);

    % 7. Use meshgrid to define a discrete grid over the state-space, and at each point of the grid calculate the matrix measure induced by the weighted Euclidean norm.
    % TODO make grid_matmis generic for any mesh grid
    %     mu_matrix_L1 = grid_matmis(A, delta_2, omega_2, 'L1');
    mu_matrix_L2 = grid_matmis(A, delta_2, omega_2, 'L2');
    MU = real(mu_matrix_L2); %remove numerical complex error 
    %     mu_matrix_Linf = grid_matmis(A, delta_2, omega_2, 'Linf');

    %find the minimal P_final-norm
    mat_P_norm = 0*delta_2;
    for row = 1:length(delta_2(:,1))
        for col = 1:length(delta_2(1,:))
            mat_P_norm(row,col) = p_norm(P_final,delta_2(row,col),omega_2(row,col),d_eq,w_eq);
        end
    end
    idx_threshold = find(MU<1e-2 & MU>-1e-2); %linear indexing
    mat_P_morm_mu0 = mat_P_norm(idx_threshold);
    min_dist = min(mat_P_morm_mu0(:));
    %{
    figure(10)
    mesh(delta_2/pi, omega_2/(2*pi),mat_P_norm);
    zlim([-5,50])
    %}
    % 8. Use contour (or any other function of your choice) to draw the boundary of the convergence, that is, the line along which the matrix measure is zero.
    figure(1)
    hold on;
    contourf(delta_2/pi, omega_2/(2*pi), -MU, [0, 0]);% minus for coloring negative values
    contour(delta_2/pi, omega_2/(2*pi), mat_P_norm, [min_dist,min_dist]);
    scatter (d_eq/pi,w_eq/(2*pi),'b*')
    axis square;
    xlabel('delta 2 [rad/\pi]')
    ylabel('omega [Hz]')
    title('Areas where \mu <0')
    

    figure (2);
    hold on;
    plot(delta_2/pi, mu_matrix_L2(1, :));
    xlabel('delta 2 [rad/\pi]')
    ylabel('matrix measure')
    title('the matrix measure of the jacovian in defferent angles')
    hold off;


%%simulate the one gen model
w_array = linspace(49.9*2*pi,50.1*2*pi,4);
d_0_array = linspace(0,0.4*pi,4);
for ij = 1:length (w_array)
    for ii = 1:length(d_0_array)
        w_0 = w_array(ij);
        d_0 = d_0_array(ii);
        [w,d2,time] = one_gen_model_sim;

        %plot all the signals in the same plot
        figure (3);
        hold on
        plot (d2/pi,w/(2*pi))
        ylabel ('omega [Hz]')
        xlabel ('delta 2 [rad/\pi]')
        %xlim ([-2,2]);
        %ylim([46,54]);

        %plot distance in P terms
        dist_w = w-w(length(w));
        dist_d2 = d2-d2(length(d2));
        norm_P = 0*w;
        for in = 1:length(w)
            norm_P(in) = p_norm(P_final,d2(in),w(in),d2(length(d2)),w(length(w)));
        end
        e_time = norm_P(1)*exp(-time);
        
        figure;
        hold on;
        plot (time,norm_P)
        plot (time,e_time)
        xlabel ('time')
        ylabel ('P norm of distance to eq. point')
        xlim ([0,10]);
        legend ('data','exponent')
        title(['start point w=',num2str(w(1)/(2*pi)),'Hz, d2=',num2str(d2(1)/pi),'rad/\pi'])
        hold off;
        %}

        %check if the exponnent is smaller than the calculated norm
        exp_smaller_than_norm = e_time<norm_P;
        figure (1)
        if (sum(exp_smaller_than_norm)>0)
            scatter(d_0/pi,w_0/(2*pi),'rx');
        else
            scatter (d_0/pi,w_0/(2*pi),'gx');
        end
    end
end
end

% for base example, don't need to solve it just return the known points
function equilibrium_points = solve_power_flow()
    global P_ref a_21 ws P_L
    equilibrium_points = cell([1, 2]);
    equilibrium_points{1} = [(asin((P_ref-P_L) / a_21)), ws]; % (delta,omega) point
    equilibrium_points{2} = [pi - (asin((P_ref-P_L) / a_21)), ws]; % (delta,omega) point
end

%{
function equilibrium_points = solve_power_flow(delta_2, omega_2)
    global P_ref a_21 P_L1 P_L2
    P_1 = -a_21 * sin(delta_2) +P_L1;
    P_2 = a_21 * sin(delta_2) +P_L2;

end
%}