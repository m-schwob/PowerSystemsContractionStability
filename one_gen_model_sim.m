clear;clc;close all;
%simulation parameters
K = 8.73e-7;
D = 0.18*K;
P_ref = 4e6;
P_L = -5e6;
ws = 50*2*pi;
a_21 = (12e3)^2 * abs(-1/(1j*ws*0.0275));

%set simulation timescale
t_final = 30;
t_max_step = 0.05;

%simulation inputs
w_array = 49:1:51;
d_0_array = -pi/2:1:pi/2;
for ij = 1:length (w_array)
for ii = 1:length(d_0_array)
w_0 = w_array(ij)*2*pi;
d_0 = d_0_array(ii);
sim_1 = sim("one_gen_model.slx");

%get data
w_signal = sim_1.yout.getElement('w');
w = w_signal.Values.data;
d2_signal = sim_1.yout.getElement('d_2');
d2 = d2_signal.Values.data;

%plot
figure (1)
hold on
plot (d2/pi,w/(2*pi))
ylabel ('omega [Hz]')
xlabel ('delta 2 [rad/\pi]')
xlim ([-2,2]);
ylim([49,51]);
end
end

%%
k_array = [1e-10,1e-9,1e-8,2e-8,7e-8,1e-7,1e-6,1e-5,5e-5,6e-5]
for ii = 1:length (k_array)
w_0 = 49*2*pi;
d_0 = -0.4*pi;
K = k_array(ii);
D = 0.18*K;
sim_1 = sim("one_gen_model.slx");

%get data
w_signal = sim_1.yout.getElement('w');
w = w_signal.Values.data;
d2_signal = sim_1.yout.getElement('d_2');
d2 = d2_signal.Values.data;

%plot
figure (2)
hold on
plot (d2/pi,w)
ylabel ('omega [rad/s]')
xlabel ('delta 2 [rad/\pi]')
xlim ([-2,2]);
ylim([49*2*pi,51*2*pi]);
end

