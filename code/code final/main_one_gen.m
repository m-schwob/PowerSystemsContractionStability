clear; clc; close all;
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

% define routs' parametes 
w_bounds = [49*2*pi,51*2*pi];
d2_bounds = [-0.1*pi,0.4*pi];
num_of_d2_points = 3;
num_of_w_points = 3;

%define Q for lyaponuv
%
Q = [1,0;0,1];
MU_lim = [-0.5e-2,0.5e-2];
color = 'r';
%{
Q = [2,-1;-1,2];
MU_lim = [-0.5e-1,0.5e-1];
color = 'c';
%
Q = [100,1;1,1];
MU_lim = [-0.5e-1,0.5e-1];
color = 'g';
%
Q = [1000,0.1;0.1,1];
MU_lim = [-0.5e-3,0.5e-3];
color = 'b';
%
Q = [100,9;9,10];
MU_lim = [-1e-2,1e-2];
color = 'm';
%}

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
    P = lyap(J_eq, Q);

    % 6. Find the square root of the result of lyap using sqrtm. The transformation we will use to define the weighted Euclidean norm is the inverse of the square root.
    T = sqrtm(P);
    T_final = T^(-1); % matlab uses a different lyapunov eq
    A = @(d) (T_final * J(d) * T_final^(-1));
    mat_P_norm = zeros(len_x,len_y);
    mu_matrix_L2 = mat_P_norm;
    temp_mat = zeros(2);
    for x = 1:len_x
        for y = 1:len_y
            mu_matrix_L2(x,y) = matmis(A(d2(y)),'L2');
            temp_mat = T_final*[w(x)-w_eq;d2(y)-d2_eq];
            mat_P_norm(x,y) = norm(temp_mat);
        end
    end
    MU = real(mu_matrix_L2); %remove numerical complex error
    idx_threshold = find(MU>MU_lim(1) & MU<MU_lim(2)); %linear indexing
    mat_P_morm_mu0 = mat_P_norm(idx_threshold); % only the places where mu=0 (approximetly)
    min_dist = min(mat_P_morm_mu0(:));

    figure(1)
    hold on;
    contourf(d2/pi,w/(2*pi),-MU, [0,0],'red');
    scatter(d2_eq/pi,w_eq/(2*pi), '*')
    contour(d2/pi, w/(2*pi), mat_P_norm, [min_dist,min_dist],'g');
    xlabel('delta 2 [rad/\pi]')
    ylabel('omega 1 [Hz]')
    title('Area where \mu <0')
    hold off;

    figure(4)
    hold on;
    scatter(d2_eq/pi,w_eq/(2*pi), 'k*')
    contour(d2/pi, w/(2*pi), mat_P_norm, [min_dist,min_dist],color);
    xlabel('delta 2 [rad/\pi]')
    ylabel('omega 1 [Hz]')
    title('Contraction Areas with Different Q Matrices')
    hold off;

    %{
    figure (2);
    hold on;
    plot(d2/pi, MU(1,:));
    xline(d2_eq/pi);
    yline(0)
    xlabel('delta 2 [rad/\pi]')
    ylabel('matrix measure')
    title('the matrix measure of the jacovian in defferent angles')
    hold off;
    %}

    %%simulate the one gen model
    w_array = linspace(w_bounds(1),w_bounds(2),num_of_w_points);
    d2_array = linspace(d2_bounds(1),d2_bounds(2),num_of_d2_points);
    for ij = 1:length (w_array)
        for ii = 1:length(d2_array)
            w_0 = w_array(ij);
            d_0 = d2_array(ii);
            [w_sim,d2_sim,time] = one_gen_model_sim;

            %plot all the signals in the same plot
            %
            figure (3);
            hold on
            plot (d2_sim/pi,w_sim/(2*pi))
            ylabel ('omega [Hz]')
            xlabel ('delta 2 [rad/\pi]')
            %xlim ([-2,2]);
            %ylim([46,54]);
            %}

            %plot distance in P terms
            norm_P = 0*w_sim;
            for in = 1:length(w_sim)
                temp_mat = T_final*[w_sim(in)-w_eq;d2_sim(in)-d2_eq];
                norm_P(in) = norm(temp_mat);
            end
            % remove all the places where norm_P<=0:
            norm_P = norm_P(norm_P>0);
            time = time(norm_P>0);
            %calculate the bound
            e_time = 0*time;
            start_dist = norm_P(1);
            for t = 1:length(time)
                idx_to_find_max = find(mat_P_norm<start_dist); %linear indexing
                etha = max(MU(idx_to_find_max));
                if (etha>0) % if etha>0 then etha doesn't mean anything
                    % and we need  to prove that the graph is not bounded
                    % by a descending exponent
                    shape = '+';
                    final_dist = min(norm_P);
                    final_time = time(end);
                    e_time(t) = (final_dist-start_dist)*time(t)/final_time + start_dist;
                else
                    shape = 'x';
                    e_time(t) = start_dist*exp(etha*time(t));
                end
            end
            %{
            figure("Name",['start point w=',num2str(w_sim(1)/(2*pi)),'Hz, d2=',num2str(d2_sim(1)/pi),'rad/\pi']);
            hold on;
            plot (time,norm_P)
            plot (time,e_time)
            xlabel ('time')
            ylabel ('P norm of distance to eq. point')
            legend ('data','exponent')
            title(['start point w=',num2str(w_sim(1)/(2*pi)),'Hz, d2=',num2str(d2_sim(1)/pi),'rad/\pi'])
            hold off
            %}

            %check if the exponnent is smaller than the calculated norm
            exp_minus_norm = e_time-norm_P;
            exp_smaller_than_norm = any(exp_minus_norm<-1e-15);
            figure(1)
            hold on

            if (shape == '+')
                scatter(d_0/pi,w_0/(2*pi),6,'r',shape);
            else
                scatter (d_0/pi,w_0/(2*pi),6,'g',shape);
            end

        end
    end
    figure(3)
    hold on
    scatter(d2_eq/pi,w_eq/(2*pi),'o')
    contour(d2/pi, w/(2*pi), mat_P_norm, [min_dist,min_dist],color);
   
end

function equilibrium_points = solve_power_flow()
global P_ref a_21 ws
equilibrium_points = cell([1, 2]);
equilibrium_points{1} = [(asin((P_ref) / a_21)), ws]; % (delta,omega) point
equilibrium_points{2} = [pi - (asin((P_ref) / a_21)), ws]; % (delta,omega) point
end