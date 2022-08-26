
%runs a simulation of the two gen model
function [w1,w2,d2,time,w1_eq,w2_eq,d2_eq] = two_gen_model_sim
    sim_1 = sim("two_gen_model_2020b.slx");
    %get data
    w1_signal = sim_1.yout{1};
    w1 = w1_signal.Values.data;
    w2_signal = sim_1.yout{2};
    w2 = w2_signal.Values.data;
    time = w1_signal.Values.time;
    d2_signal = sim_1.yout{3};
    d2 = d2_signal.Values.data;
    w1_eq = w1(end);
    w2_eq = w2(end);
    d2_eq = d2(end);
end
