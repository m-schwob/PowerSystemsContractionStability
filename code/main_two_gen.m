clear; clc;
%  SYSTEM DEFINITIONS  %

% define system constants
global K1 D1 P1_ref K2 D2 P2_ref a_21 ws
K1 = 8.73e-8; %0.7 works best
D1 = 0.18 * K1;
P1_ref = 1e7;
K2 = 6.12e-8; %0.7 works best
D2 = 0.18 * K1;
P2_ref = 1e7;
ws = 50 * 2 * pi;
a_21 = (12e3)^2 * abs(1 / (1j * ws * 0.0275));



% define boundaries
w1_min = 49 * 2 * pi; % rad/sec
w1_max = 51 * 2 * pi; % rad/sec
w2_min = 49 * 2 * pi; % rad/sec
w2_max = 51 * 2 * pi; % rad/sec
d2_min = -pi / 2; % rad
d2_max = pi / 2; % rad

%simulation parameters
t_final = 10;
t_max_step = 0.01;

% define resolutions
w1_res = 0.1;
w2_res = 0.1;
d2_res = 0.1;
w1_bounds = [49.7*2*pi,50.3*2*pi];
w2_bounds = [49.7*2*pi,50.3*2*pi];
d2_bounds = [-0.4*pi,0.4*pi];
num_of_d2_points = 6;
num_of_w1_points = 6;
num_of_w2_points = 6;

% define Jacobian formula
J = @(d2) [-K1/D1, 0, 3 * K1 * a_21 * cos(d2);
    0, -K2/D2, -3 * K1 * a_21 * cos(d2);
    -1,1,0]; % top to bottom: (w1,w2,d2)

% define numerical grid
d2 = d2_min:d2_res:d2_max;
w1 = w1_min:w1_res:w1_max;
w2 = w2_min:w2_res:w2_max;
len_x = length(w1);
len_y = length(w2);
len_z = length(d2);
%  ALGORITHM  %

% run one sim to find the equalibriom point
w1_array = linspace(w1_bounds(1),w1_bounds(2),num_of_w1_points);
w2_array = linspace(w2_bounds(1),w2_bounds(2),num_of_w2_points);
d2_array = linspace(d2_bounds(1),d2_bounds(2),num_of_d2_points);

w1_0 = w1_array(1);
w2_0 = w2_array(1);
d2_0 = d2_array(1);
[~,~,~,~,w1_eq,w2_eq,d2_eq] = two_gen_model_sim;
J_eq = J(d2_eq);

%calculate the Mu matrix using liaponov
Q = eye(3);
P = lyap(J_eq, Q);
T = sqrtm(P);
T_final = T^-1; % matlab uses a different lyapunov eq
A = @(d) (T_final * J(d) * T_final^(-1));
mat_P_norm = zeros(len_x,len_y,len_z);
mu_matrix_L2 = mat_P_norm;
temp_mat = zeros(3);
for x = 1:len_x
    for y = 1:len_y
        for z = 1:len_z
            mu_matrix_L2(x,y,z) = matmis(A(d2(z)),'L2');
            temp_mat = T_final*[w1(x)-w1_eq;w2(y)-w2_eq;d2(z)-d2_eq];
            mat_P_norm(x,y,z) = norm(temp_mat);
        end
    end
end
MU = real(mu_matrix_L2); %remove numerical complex error
MU_2d = squeeze(MU(:,1,:)); % only w1 and d2
mat_P_norm_2d = squeeze(mat_P_norm(:,1,:));
idx_threshold = find(MU<1e-1 & MU>-1e-1); %linear indexing
mat_P_morm_mu0 = mat_P_norm(idx_threshold);
min_dist = min(mat_P_morm_mu0(:));
idx_to_dist = find(mat_P_norm<=min_dist);
[idx_x,idx_y,idx_z] = ind2sub([len_x,len_y,len_z],idx_to_dist);
w1_dist = w1(idx_x);
w2_dist = w2(idx_y);
d2_dist = d2(idx_z);
a = linspace(0,100,length(d2_dist));
%{
figure(1)
hold on;
contourf(d2/pi,w1/(2*pi),-MU_2d, [0,0],'red');
contour(d2/pi, w1/(2*pi), mat_P_norm_2d, [min_dist,min_dist],'red')
scatter(d2_eq/pi,w1_eq/(2*pi),'*')
xlabel('delta 2 [rad/\pi]')
ylabel('omega 1 [Hz]')
title('Areas where \mu <0')
hold off;
%}

%
figure (2);
hold on;
scatter3(w1_dist/(2*pi),w2_dist/(2*pi),d2_dist/pi,40,a,"o","filled",'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',0)
zlabel('delta 2 [rad/\pi]')
ylabel('omega 2 [Hz]')
xlabel('omega 1 [Hz]')
title({'Areas where the distance in T term', 'is less than the minimal distance to the line where \mu=0'})
%}

% run the rest of the simulations
for ij = 1:length(w1_array)
    for ii = 1:length(d2_array)
        for ik = 1:length(w2_array)
            w1_0 = w1_array(ij);
            w2_0 = w2_array(ik);
            d2_0 = d2_array(ii);
            dist_matrix = T_final*[w1_0-w1_eq;w2_0-w2_eq;d2_0-d2_eq];
            P_norm_0 = norm(dist_matrix);
            if (P_norm_0<min_dist)
                in_area = true;
            else
                in_area = false;
            end
            [w1_sim,w2_sim,d2_sim,time,w1_eq_sim,w2_eq_sim,d2_eq_sim] = two_gen_model_sim;

            %plot all the signals in the same plot
            %
            figure (3);
            hold on
            grid on
            plot3 (w1_sim/(2*pi),w2_sim/(2*pi), d2_sim/pi)
            xlabel ('omega 1 [Hz]')
            ylabel ('omega 2 [Hz]')
            zlabel ('delta 2 [rad/\pi]')
            %}

            %plot distance in P terms
            norm_P = 0*w1_sim;
            for in = 1:length(w1_sim)
                temp_mat = T_final*[w1_sim(in)-w1_eq;w2_sim(in)-w2_eq;d2_sim(in)-d2_eq];
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
                    final_dist = min(norm_P);
                    final_time = time(end);
                    e_time(t) = (final_dist-start_dist)*time(t)/final_time + start_dist;
                else
                    e_time(t) = start_dist*exp(etha*time(t));
                end
            end
            %{
            figure;
            hold on;
            plot (time,norm_P)
            plot (time,e_time)
            xlabel ('time')
            ylabel ('P norm of distance to eq. point')
            xlim ([0,10]);
            legend ('data','exponent')
            title(['start point w1=',num2str(w1_sim(1)/(2*pi)),'Hz, w2=',num2str(w2_sim(1)/(2*pi)),'Hz, d2=',num2str(d2_sim(1)/pi),'rad/\pi'])
            hold off
            %}

            %check if the exponnent is smaller than the calculated norm
            exp_minus_norm = e_time-norm_P;
            exp_smaller_than_norm = any(exp_minus_norm<-1e-3);
            if (exp_smaller_than_norm)
                fprintf("exp smaller than norm")
            end
            figure(2)
            hold on
            if(in_area)
                scatter3(w1_0/(2*pi),w2_0/(2*pi),d2_0/pi,30,'green','filled','o','MarkerEdgeColor','k');
            else
                scatter3(w1_0/(2*pi),w2_0/(2*pi),d2_0/pi,30,'red','*');
            end
        end
    end
end
figure(3)
scatter3(w1_eq/(2*pi),w2_eq/(2*pi),d2_eq/pi,'green','filled','o');
figure(2)
scatter3(w1_eq/(2*pi),w2_eq/(2*pi),d2_eq/pi,40,'green','filled','diamond','MarkerEdgeColor','k');

