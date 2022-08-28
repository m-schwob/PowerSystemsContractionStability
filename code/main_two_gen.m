clear; close all; clc;
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
w1_array = linspace(49.9*2*pi,50.1*2*pi,3);
w2_array = linspace(49.9*2*pi,50.1*2*pi,3);
d2_0_array = linspace(0,0.4*pi,3);

w1_0 = w1_array(1);
w2_0 = w2_array(1);
d2_0 = d2_0_array(1);
[~,~,~,~,w1_eq,w2_eq,d2_eq] = two_gen_model_sim;
J_eq = J(d2_eq);

%calculate the Mu matrix using liaponov
Q = eye(3);
P = lyap(J_eq, Q);
T = sqrtm(P);
P_final = T^-1;
A = @(d) (T^(-1) * J(d) * T);
mat_P_norm = zeros(len_x,len_y,len_z);
mu_matrix_L2 = mat_P_norm;
temp_mat = zeros(3);
for x = 1:len_x
    for y = 1:len_y
        for z = 1:len_z
            mu_matrix_L2(x,y,z) = matmis(A(d2(z)),'L2');
            temp_mat = P_final*[w1(x)-w1_eq;w2(y)-w2_eq;d2(z)-d2_eq];
            mat_P_norm(x,y,z) = norm(temp_mat);
        end
    end
end
MU = real(mu_matrix_L2); %remove numerical complex error
MU_2d = squeeze(MU(1,:,:));
idx_threshold = find(MU<1e-1 & MU>-1e-1); %linear indexing
mat_P_morm_mu0 = mat_P_norm(idx_threshold);
min_dist = min(mat_P_morm_mu0(:));
idx_to_dist = find(mat_P_norm<=min_dist);
[idx_x,idx_y,idx_z] = ind2sub([len_x,len_y,len_z],idx_to_dist);
w1_dist = w1(idx_x);
w2_dist = w2(idx_y);
[~,d2_dist] = meshgrid(d2(idx_z),d2(idx_z));

figure(1)
hold on;
contourf(d2/pi,w1/(2*pi),-MU_2d, [0,0],'red');
scatter(d2_eq/pi,w1_eq/(2*pi),'*')
xlabel('delta 2 [rad/\pi]')
ylabel('omega 1 [Hz]')
title('Areas where \mu <0')
hold off;

%{
figure (2);
hold on;
surf(w1_dist/(2*pi),w2_dist/(2*pi),d2_dist/pi);
%}

% run the rest of the simulations
for ij = 1:length (w1_array)
    for ii = 1:length(d2_0_array)
        for ik = 1:length(w2_array)
            w1_0 = w1_array(ij);
            w2_0 = w2_array(ik);
            d2_0 = d2_0_array(ii);
            [w1_sim,w2_sim,d2_sim,time,w1_eq_sim,w2_eq_sim,d2_eq_sim] = two_gen_model_sim;
            % TODO chek eq is eq

            %plot all the signals in the same plot
            figure (3);
            hold on
            grid on
            plot3 (w1_sim/(2*pi),w2_sim/(2*pi), d2_sim/pi)
            xlabel ('omega 1 [Hz]')
            ylabel ('omega 2 [Hz]')
            zlabel ('delta 2 [rad/\pi]')
            %xlim ([-2,2]);
            %ylim([46,54]);
            %plot distance in P terms

            norm_P = 0*w1_sim;
            for in = 1:length(w1_sim)
                temp_mat = P_final*[w1_sim(in)-w1_eq;w2_sim(in)-w2_eq;d2_sim(in)-d2_eq];
                norm_P(in) = norm(temp_mat);
            end
            %calculate the bound
            e_time = 0*time;
            for t = 1:length(time)
                integral = 0;
                for t_in = 1:t
                    [~,idx_d2] = min(abs(d2_sim(t_in)-d2));
                    [~,idx_w1] = min(abs(w1_sim(t_in)-w1));
                    [~,idx_w2] = min(abs(w2_sim(t_in)-w2));
                    if (t_in ~=1)
                        integral = integral + (time(t_in)-time(t_in-1))*200*(MU(idx_w1,idx_w2,idx_d2)-MU(idx_w1_prev,idx_w2_prev,idx_d2_prev));
                    end
                    idx_w1_prev = idx_w1;
                    idx_w2_prev = idx_w2;
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
            title(['start point w1=',num2str(w1_sim(1)/(2*pi)),'Hz, w2=',num2str(w2_sim(1)/(2*pi)),'Hz, d2=',num2str(d2_sim(1)/pi),'rad/\pi'])
            hold off
            %}

            %check if the exponnent is smaller than the calculated norm
            exp_smaller_than_norm = e_time<norm_P;
            figure(1)
            hold on
            if (sum(exp_smaller_than_norm)>0)
                scatter(d2_0/pi,w1_0/(2*pi),'rx');
            else
                scatter (d2_0/pi,w1_0/(2*pi),'gx');
            end
        end
    end
end
figure(3)
scatter3(w1_eq/(2*pi),w2_eq/(2*pi),d2_eq/pi,'ro');

