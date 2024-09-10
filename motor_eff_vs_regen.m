%% Motor Efficiency VS Efficiency Points
clear; clc

options.t_driverchange = 180;
options.fastest_laptime = 91;
options.total_laps = 16;
options.dist_8laps = 11;
options.dist_16laps = 22;
options.lap_dist = 1.375;                       % 22 km in 16 laps = 1.375 km per lap
options.end_time_8laps = 1e4;
options.regen_favorable = false;
options.min_endur_t_specific_cell = 1.4650e3;   % only for MY21pack()... for other cells, see lines 111-112 in BatteryPointsModel
battery = MY21pack();

P_draw_cmd_no_regen =  25e3;        % no regen power [Watts]
P_draw_cmd_regen = 0;               % regen power [kWh]
P_draw_cmd_regen = P_draw_cmd_regen / ( battery.J_to_kWh * (options.total_laps * options.fastest_laptime) ); % regen power [Watts]

sim = BatterySimClass();
[time, state, log] = sim.run_endurance_thermals(P_draw_cmd_no_regen, P_draw_cmd_regen, battery, options);

event = 'efficiency';
E_regen = P_draw_cmd_regen * (options.total_laps * options.fastest_laptime) * battery.J_to_kWh;

ratio = 1 ./ efficiency;
% EXPLANATION OF EYOUR EQN BELOW: (https://wikis.mit.edu/confluence/display/FSAECONTROLS/Motor+Efficiency+VS+Battery+Efficiency+Points)
% tmin = fastest endurance time recorded at competition (out of all competing cars) [seconds]
% Emin = least amount of net energy used at competition (out of all competing cars) [kWh]
% tyour = endurance time (our car) [seconds]
% Eyour = net energy used (our car) [kWh]
tmin = options.min_endur_t_specific_cell;
Emin = 1.9; % Real minimum energy use based on MY19
tyour = ones(1,length(state))*time(end);
Eyour = ( (state(end,3)*battery.J_to_kWh + E_regen) .* (ratio) ) - ( E_regen ./ (ratio.^2) );

efficiency = linspace(0,1, length(state));
pts = dynamic_pts(event, tmin, tyour, Emin, Eyour);

figure;
plot(efficiency, pts)
xlabel('Efficiency')
ylabel('Points')
xlim([.7,1])
improvePlot

figure;
plot(efficiency, Eyour)
xlabel('Efficiency')
ylabel('Net Energy (kWh)')
xlim([.7,1])
improvePlot