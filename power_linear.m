P = (1:30)*1000;
t = find_laptime_from_power(P);
d = 22000/16;
v = d./t;

figure;plot(P,v)
hold on;plot([P(1) P(end)], [v(1) v(end)])
improvePlot
%%

for i=1:length(P)
    P_avg = P(i);
    P_rms = P(i) + 28000;
    t_correct = 22000 / get_v(P_rms,P_avg);
    t_sim =     22000 / get_v(P_avg,P_avg);
    delta(i) = abs(t_correct - t_sim);
end

max(delta)

function v = get_v(P_rms, P_avg)
    a = 8.418e-10;
    b = 8.03e-5;
    c = 11.94;
    v = a*P_rms + b*P_avg + c;
end