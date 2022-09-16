
%runs a simulation of the one gen model
function [w,d2,time] = one_gen_model_sim
    sim_1 = sim("one_gen_model_2020b.slx");

    %get data
    w_signal = sim_1.yout{1};
    w = w_signal.Values.data;
    time = w_signal.Values.time;
    d2_signal = sim_1.yout{2};
    d2 = d2_signal.Values.data;
end
