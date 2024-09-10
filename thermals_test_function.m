clear;clc

sim = BatterySimClass();
options.t_driverchange = 180;
options.fastest_laptime = 91;
options.total_laps = 16;
options.dist_8laps = 11;
options.dist_16laps = 22;
options.lap_dist = 1.375;   % 22 km in 16 laps = 1.375 km per lap
options.end_time_8laps = 1e4;
options.regen_favorable = false;
options.min_endur_t_specific_cell = 1.4650e3; % only for MY21pack()... for other cells, see lines 111-112 in BatteryPointsModel
battery = MY21pack();

P_draw_no_regen_vec =  [25e3] ; % non regen power (Watts)
P_draw_cmd_regen = 0 / (battery.J_to_kWh * (options.total_laps * options.fastest_laptime) ) ; % regen power (Watts)

for i = 1:length(P_draw_no_regen_vec)
    
    [time, state, log] = sim.run_endurance_thermals(P_draw_no_regen_vec(i), P_draw_cmd_regen, battery, options);

    %close all
    figure
    subplot(121)
    hold all
    plot(time, state(:,1),'k','DisplayName','Temp (C)')
    state(:,1)
    plot(time, [log.voltage]/10, 'c','DisplayName','Voltage (V/10)')
    plot(time, [log.ocv]/10,'c--','DisplayName','OCV (V/10)')

    plot(time, [log.P_draw_no_regen]/1e3,'r','DisplayName','P draw accel (kW)')
    plot(time, [log.P_draw_cmd_no_regen]/1e3,'r--','DisplayName','P draw cmd accel (kW)')
    plot(time, [log.P_draw_regen]/1e3,'b','DisplayName','P draw regen (kW)')
    plot(time, [log.P_draw_cmd_regen]/1e3,'b--','DisplayName','P draw cmd regen (kW)')
    ylim([0 100])
    legend

    subplot(122)
    hold all
    plot(time, state(:,2),'g','DisplayName','Distance (km)')
    plot(time, [log.laptime]/10,'k', 'DisplayName','Laptime (s/10)')
    text(time(end),state(end,2), [num2str(state(end,2)) ' km'])

    plot(time, state(:,3)*battery.J_to_kWh,'c','DisplayName','Energy Used (kWh)')
    text(time(end),state(end,3)*battery.J_to_kWh, [num2str(state(end,3)/36e5) ' kWh'])

    plot(time, [log.P_diss]/1e3, 'r','DisplayName','P diss (kW)')
    plot(time, [log.P_cool]/1e3, 'b','DisplayName','P cool (kW)')

    legend
    improvePlot

    figure
    hold all
    plot(time, [log.I_draw_regen], 'DisplayName', 'Regen Current')
    plot(time, [log.I_draw_no_regen], 'DisplayName', 'NonRegen Current')
    I_draw = [log.I_draw_regen] + [log.I_draw_no_regen];
    plot(time, I_draw, 'DisplayName', 'Total Current Draw')

    xlabel('Time')
    ylabel('Current')
    legend
    improvePlot
    
end