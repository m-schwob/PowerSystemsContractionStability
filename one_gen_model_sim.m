clear;clc;close all;
%simulation parameters
K = 8.73e-7;
D = 0.18*K;
P_ref = 4e6;
P_L = P_ref;
ws = 50*2*pi;
a_21 = (12e3)^2 * abs(-1/(1j*ws*0.0275));

%set simulation timescale
t_final = 10;
t_max_step = 0.05;

%simulation inputs
w_0 = 48*2*pi;
d_0 = pi/4;
sim_1 = sim("one_gen_model.slx");

%get data
w_signal = sim_1.yout.getElement('w');
w = w_signal.Values.data;
d2_signal = sim_1.yout.getElement('d_2');
d2 = d2_signal.Values.data;

figure (1)
hold on
plot (d2,w)
ylabel ('omega [rad/s]')
xlabel ('delta 2 [rad]')



