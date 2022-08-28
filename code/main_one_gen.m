clear; close all; clc;
%  SYSTEM DEFINITIONS  %

% define system constants
global K D a_21 P_ref ws
K = 8.73e-8; %0.7 works best
D = 0.18 * K;
P_ref = 1e7;
ws = 50 * 2 * pi;
a_21 = (12e3)^2 * abs(1 / (1j * ws * 0.0275));



% define boundaries
w_min = 49 * 2 * pi; % rad/sec
w_max = 51 * 2 * pi; % rad/sec
d2_min = -pi / 2; % rad
d2_max = pi / 2; % rad

%simulation parameters
t_final = 10;
t_max_step = 0.01;

% define resolutions
w_res = 0.01;
d2_res = 0.01;

% define Jacobian formula
J = @(d) [ -K / D, -3 * K * a_21 * cos(d);1 0]; % (w,d)

% define numerical grid
d2 = d2_min:d2_res:d2_max;
w = w_min:w_res:w_max;
len_x = length(w);
len_y = length(d2);
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
    d2_eq = equilibrium_points{i}(1);
    w_eq = equilibrium_points{i}(2);
    Q = eye(ndims(J_eq));
    P = lyap(J_eq, Q);

    % 6. Find the square root of the result of lyap using sqrtm. The transformation we will use to define the weighted Euclidean norm is the inverse of the square root.
    T = sqrtm(P);
    P_final = T^-1;
    A = @(d) (T^(-1) * J(d) * T);
    mat_P_norm = zeros(len_x,len_y);
    mu_matrix_L2 = mat_P_norm;
    temp_mat = zeros(2);
    for x = 1:len_x
        for y = 1:len_y
            mu_matrix_L2(x,y) = matmis(A(d2(y)),'L2');
            temp_mat = P_final*[w(x)-w_eq;d2(y)-d2_eq];
            mat_P_norm(x,y) = norm(temp_mat);
        end
    end
    MU = real(mu_matrix_L2); %remove numerical complex error
    idx_threshold = find(MU<1e-2 & MU>-1e-2); %linear indexing
    mat_P_morm_mu0 = mat_P_norm(idx_threshold);
    min_dist = min(mat_P_morm_mu0(:));

    figure(1)
    hold on;
    contourf(d2/pi,w/(2*pi),-MU, [0,0],'red');
    scatter(d2_eq/pi,w_eq/(2*pi), '*')
    contour(d2/pi, w/(2*pi), mat_P_norm, [min_dist,min_dist],'red');
    xlabel('delta 2 [rad/\pi]')
    ylabel('omega 1 [Hz]')
    title('Area where \mu <0')
    hold off;

    figure (2);
    hold on;
    plot(d2/pi, MU(1,:));
    xline(d2_eq/pi);
    yline(0)
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
            [w_sim,d2_sim,time] = one_gen_model_sim;

            %plot all the signals in the same plot
            figure (3);
            hold on
            plot (d2_sim/pi,w_sim/(2*pi))
            ylabel ('omega [Hz]')
            xlabel ('delta 2 [rad/\pi]')
            %xlim ([-2,2]);
            %ylim([46,54]);

            %plot distance in P terms

            norm_P = 0*w_sim;
            for in = 1:length(w_sim)
                temp_mat = P_final*[w_sim(in)-w_eq;d2_sim(in)-d2_eq];
                norm_P(in) = norm(temp_mat);
            end
            %calculate the bound
            e_time = 0*time;
            for t = 1:length(time)
                integral = 0;
                for t_in = 1:t
                    [~,idx_d2] = min(abs(d2_sim(t_in)-d2));
                    [~,idx_w] = min(abs(w_sim(t_in)-w));
                    if (t_in ~=1)
                    [~,idx] = min(norm_P(norm_P>0));    
                    integral = integral + (time(t_in)-time(t_in-1))*200*(MU(idx_w,idx_d2)-MU(idx_w_prev,idx_d2_prev));
                    end
                    idx_w_prev = idx_w;
                    idx_d2_prev = idx_d2;
                end
                e_time(t) = norm_P(1)*exp(integral*time(t));
            end
            %
            figure;
            hold on;
            plot (time,norm_P)
            plot (time,e_time)
            xlabel ('time')
            ylabel ('P norm of distance to eq. point')
            xlim ([0,10]);
            legend ('data','exponent')
            title(['start point w=',num2str(w_sim(1)/(2*pi)),'Hz, d2=',num2str(d2_sim(1)/pi),'rad/\pi'])
            hold off
            %}

            %check if the exponnent is smaller than the calculated norm
            exp_smaller_than_norm = e_time<norm_P;
            figure(1)
            hold on
            
            if (sum(exp_smaller_than_norm)>0)
                scatter(d_0/pi,w_0/(2*pi),'rx');
            else
                scatter (d_0/pi,w_0/(2*pi),'gx');
            end
        end
    end
    figure(3)
    hold on
    scatter(d2_eq/pi,w_eq/(2*pi),'o')
end

function equilibrium_points = solve_power_flow()
global P_ref a_21 ws
equilibrium_points = cell([1, 2]);
equilibrium_points{1} = [(asin((P_ref) / a_21)), ws]; % (delta,omega) point
equilibrium_points{2} = [pi - (asin((P_ref) / a_21)), ws]; % (delta,omega) point
end