function plot_graph_dist_two_gen(w1_0,w2_0,d2_0,T_final,w1_eq,w2_eq,d2_eq,mat_P_norm,MU)
[w1_sim,w2_sim,d2_sim,time,~,~,~] = two_gen_model_sim;

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
    etha = max(MU(mat_P_norm<start_dist));
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

figure(8);
hold on;
plot (time,norm_P)
plot (time,e_time)
xlabel ('time')
ylabel ('P norm of distance to eq. point')
title(['start point w1=',num2str(w1_sim(1)/(2*pi)),'Hz, w2=',num2str(w2_sim(1)/(2*pi)),'Hz, d2=',num2str(d2_sim(1)/pi),'rad/\pi'])
hold off

end