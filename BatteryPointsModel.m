clear
close all
clc;

%% Control which sensitivity you want to sweep, regen or power efficiency
regen_sens = true;
eff_sens = ~regen_sens;

%% Choose which pack parameters you want

%battery = MY18pack();
%battery = MY20pack();
battery = MY21pack();
%battery = GermanPouch();
%battery = samsung30q();
%battery = samsung30T();
%battery = samsung40T();
%battery = SonyVTC6();
%battery = SonyVTC6A();

%% Setup constants and input values

options.t_driverchange = 180; % [seconds]
options.fastest_laptime = 91;
options.total_laps = 16;
options.dist_8laps = 11; % [kilometers]
options.dist_16laps = 22; % [kilometers]
options.lap_dist = 1.375;   % 22 km in 16 laps = 1.375 [km per lap]
options.end_time_8laps = 1e4;
options.min_endur_t_all_cells = 1.3695e3;
options.min_endur_t_specific_cell = 9e10;
options.regen_favorable = false;
options.regen_eff_sens = 'eff_sens';
% used for eff_sens
if strcmp(options.regen_eff_sens, 'eff_sens')
    min_eff = 0.75;
    max_eff = 0.93;
    options.eff_vec = min_eff:0.02:max_eff;
end

% Sweeping values
sweep_size = 10;
min_p_count = 4;
max_p_count = 12;
p_counts = min_p_count : max_p_count;
num_pcounts = length(p_counts);

P_regen = linspace(0,90,sweep_size); %Make sure this starts from 0
cont_power_vec = getAvgPw(P_regen)*1000; %Convert kW to W
regen_E_vec = kwhRecovered(P_regen);
no_regen_power = cont_power_vec(1);
Emin = 1.9*1/battery.J_to_kWh; %Real minimum energy use based on MY19, Joules

colors_mat = colormap(jet(length(P_regen)));
colors = mat2cell(colors_mat, ones(1,sweep_size), 3); %Look at matlab docs if you're confused

z = zeros(sweep_size,sweep_size,num_pcounts);
[points_matrix, mass_points_matrix, power_points_matrix, eff_points_matrix, ...
    t_to_60_matrix, p_count_matrix, endurance_time_matrix, unlimped_laps_matrix] = deal(z);

%% Best cell for each battery capacity
% ONLY RUN THIS SECTION IF YOU WANT TO COMPARE DIFFERENT CELLS. OTHERWISE, COMMENT OUT AND THEN RUN THE SCRIPT

% cells = {samsung30q(), samsung30T(), samsung40T(), SonyVTC6(), SonyVTC6A(), MY21pack()};
% capacity_vec = ones(1,num_pcounts);
% max_points = ones(1,num_pcounts);
% 
% % Find minimum laptime (across all cells) - used as a reference point for cell comparison
% % Use default value if cells are: {samsung30q(), samsung30T(), samsung40T(), SonyVTC6(), SonyVTC6A(), MY21pack()}
% % Otherwise, uncomment and run the next 2 lines:
% % min_endur_t = find_min_laptime_all_cells(num_pcounts, p_counts, regen_E_vec, min_p_count, cont_power_vec, no_regen_power, sweep_size, kWh_2_joules, options);
% % options.min_endur_t_all_cells = min_endur_t;
% 
% for j = 1:length(cells)
%     battery = cells{j};
%     for p = 1:num_pcounts
%         p_count = p_counts(p);
%         for i=1:length(regen_E_vec)
%             eff = 0.9;
%             [endur_times, E_used_total, P_draw_no_regen_vec, P_draw_limped, time_to_60C, mass_delta, unlimped_laps] = find_point_parameters_new(p_count, min_p_count, cont_power_vec(i), no_regen_power, regen_E_vec(i)*kWh_2_joules, i, eff, battery, sweep_size, options);
%             [points, mass_points, limp_power_points, eff_points] = delta_pts_battery(mass_delta, endur_times, options.min_endur_t_all_cells, E_used_total/kWh_2_joules, Emin/kWh_2_joules);
%             
%             points_matrix(i,:,p)         = points;
%             mass_points_matrix(i,:,p)    = mass_points;
%             power_points_matrix(i,:,p)   = limp_power_points;
%             eff_points_matrix(i,:,p)     = eff_points;
%             t_to_60_matrix(i,:,p)        = time_to_60C;
%             p_count_matrix(i,:,p)        = p_count;
%             endurance_time_matrix(i,:,p) = endur_times;
%             unlimped_laps_matrix(i,:,p)  = unlimped_laps;
%         end
%         max_points(p) = max(points_matrix(:,:,p), [], 'all');
%         capacity_vec(p) = battery.capacity;
%     end
%     plot(capacity_vec, max_points)
%     hold on
% end
% 
% xlabel('Battery Capacity')
% ylabel('Points')
% legend('samsung30q', 'samsung30T', 'samsung40T', 'SonyVTC6', 'SonyVTC6A', 'MY21pack')
% improvePlot

%% Model points, etc

% Find minimum laptime (for specific cell)
% If cell is MY21pack() and regen_sens, uncomment the next line:
% options.min_endur_t_specific_cell = 1.4650e3;
% Otherwise, run the next 2 lines:
min_endur_t = find_min_laptime_specific_cell(num_pcounts, p_counts, regen_E_vec, min_p_count, cont_power_vec, no_regen_power, sweep_size, 1/battery.J_to_kWh, battery, options);
options.min_endur_t_specific_cell = min_endur_t;

for p = 1:num_pcounts
    p_count = p_counts(p);
    if regen_sens
        for i=1:length(regen_E_vec)
            eff = 0.9;
            [endur_times, E_used_total, P_draw_no_regen_vec, P_draw_limped, time_to_60C, mass_delta, unlimped_laps] = find_point_parameters_new(p_count, min_p_count, cont_power_vec(i), no_regen_power, regen_E_vec(i)*1/battery.J_to_kWh, i, eff, battery, sweep_size, options);
            [points, mass_points, limp_power_points, eff_points] = delta_pts_battery(mass_delta, endur_times, options.min_endur_t_specific_cell, E_used_total * battery.J_to_kWh, Emin * battery.J_to_kWh);
            
            points_matrix(i,:,p)         = points;
            mass_points_matrix(i,:,p)    = mass_points;
            power_points_matrix(i,:,p)   = limp_power_points;
            eff_points_matrix(i,:,p)     = eff_points;
            t_to_60_matrix(i,:,p)        = time_to_60C;
            p_count_matrix(i,:,p)        = p_count;
            endurance_time_matrix(i,:,p) = endur_times;
            unlimped_laps_matrix(i,:,p)  = unlimped_laps;
        end
    elseif eff_sens
        for i=1:length(eff_vec)
            eff = eff_vec(i);
            [endur_times, E_used_total, P_draw, ~, time_to_60C, mass_delta] = find_point_parameters_new(p_count, min_p_count, cont_power_vec(1), no_regen_power, regen_E_vec(1)*1/battery.J_to_kWh, 1, eff, battery, sweep_size, options);
                                                                                                                                                    %i                               i                i
            [points, mass_points, limp_power_points, eff_points] = delta_pts_battery(mass_delta, endur_times, options.min_endur_t_specific_cell, E_used_total * battery.J_to_kWh, Emin * battery.J_to_kWh);
            points_matrix(i,:,p) = points;
            mass_points_matrix(i,:,p) = mass_points;
            power_points_matrix(i,:,p) = limp_power_points;
            eff_points_matrix(i,:,p) = eff_points;
            t_to_60_matrix(i,:,p) = time_to_60C;
            p_count_matrix(i,:,p) = p_count;
        end
    end
end
%% Get best parameters across all pcounts
[points_matrix_best, points_index] = max(points_matrix,[],3, 'linear'); % best points across all pcounts for each no_regen power and regen power value
mass_points_matrix_best    = mass_points_matrix(points_index);
power_points_matrix_best   = power_points_matrix(points_index);
eff_points_matrix_best     = eff_points_matrix(points_index);
t_to_60_matrix_best        = t_to_60_matrix(points_index);
p_count_matrix_best        = p_count_matrix(points_index);
endurance_time_matrix_best = endurance_time_matrix(points_index);

%% Make contour plots
x = P_draw_no_regen_vec/10^3;
x_name = 'Average Power: Acceleration (kW)';

if regen_sens
    y = P_regen;
    y_name = 'Peak Power: Regen (kW)';
elseif eff_sens
    y = eff_vec;
    y_name = 'Efficiency (%)';
end

z        = points_matrix_best;
z(:,:,2) = mass_points_matrix_best;
z(:,:,3) = power_points_matrix_best;
z(:,:,4) = eff_points_matrix_best;
z(:,:,5) = t_to_60_matrix_best;
z(:,:,6) = p_count_matrix_best;
z(:,:,7) = endurance_time_matrix_best;
z(:,:,8) = unlimped_laps_matrix(points_index);


z_names = { 'Overall Points',       ...
            'Mass Points',          ...
            'Endurance Points',     ...
            'Efficiency Points',    ...
            'Time to 60C',          ...
            'P Count',              ...
            'Endurance Time',       ...
            'Unlimped Laps'};

interpolations = ones(1,size(z,3))* 75/sweep_size;
interpolations(6) = 1; % don't interpolate the pcount graph
plot_title = ['p-count ', num2str(min_p_count),'-',num2str(max_p_count)];

multiple_double_sens_plots(x, x_name, y, y_name, z, z_names, interpolations, plot_title)

%% Get best points for each pcount to compare pcounts
[points_best_ea_pcount, points_index_ea_pcount] = max(points_matrix,[],[1 2], 'linear'); % best points across all no_regen and regen power values, for each pcount
[y_idxs,x_idxs,~] = ind2sub(size(points_matrix),points_index_ea_pcount);

regen_best_ea_pcount = y(y_idxs(:));
no_regen_best_ea_pcount = x(x_idxs(:));

points_best_ea_pcount = points_best_ea_pcount(:)-min(points_best_ea_pcount);
[best_pts, best_pts_idx] = max(points_best_ea_pcount);
best_pcount = p_counts(best_pts_idx);

figure
hold all
plot(p_counts, points_best_ea_pcount,'DisplayName','Points - Higher is better')
plot(p_counts, regen_best_ea_pcount,'DisplayName','Regen Power (kW)')
plot(p_counts, no_regen_best_ea_pcount,'DisplayName','Accel (No Regen) Power (kW)')

plot(best_pcount,best_pts,'bp','MarkerSize',20,'MarkerFaceColor','b','HandleVisibility','off')
text(best_pcount-.5,best_pts-5,['Optimal pcount ' num2str(best_pcount)])

ylim([0 100])
xlim([min(p_counts)-1,max(p_counts)+1])
xlabel('Pcount')
ylabel('Optimal for each pcount')
legend
improvePlot

%% Best Points vs P_draw_no_regen (for each p_count)
figure
hold all

temp_matrix = ones(1, sweep_size, length(p_counts));

for i = 1:length(p_counts)
    points_matrix_per_pcount = points_matrix(:,:,i);
    temp_matrix(1, :, i) = max(points_matrix_per_pcount, [], 1);
end

pos_temp_matrix = temp_matrix - min(temp_matrix, [], 'all');

for i = 1:length(p_counts)
    plot(x, pos_temp_matrix(1,:,i), 'DisplayName', ['pcount = ' num2str(p_counts(i)) ])
end

xlabel('Average Power: Acceleration (kW)')
ylabel('Points')
legend
improvePlot

%%

figure
hold all
index = 5;

a = max(points_matrix(:,index:sweep_size,1), [], 'all');
c = max(power_points_matrix(:,index:sweep_size,1), [], 'all');
d = max(eff_points_matrix(:,index:sweep_size,1), [], 'all');

maxi = max([a, c, d]);

pos_temp_matrix = points_matrix - maxi;
% colors = ['y', 'm', 'c', 'r', 'g', 'b', 'w', 'k'];
colors = {'#4DBEEE', '#D95319', '#EDB120', '#77AC30', '#7E2F8E', '#0072BD'};

for i = index:sweep_size
    plot(regen_E_vec, pos_temp_matrix(:,i,1), 'Color', colors{i-index+1}, 'DisplayName', ['Avg Accel Pwr = ' num2str(P_draw_no_regen_vec(i)/1000) ])
%     plot(regen_E_vec, mass_points_matrix(:,i,1) - maxi, '--', 'Color', colors{i-index+1}, 'Marker', 'o', 'DisplayName', 'Mass Points')
    plot(regen_E_vec, power_points_matrix(:,i,1) - maxi, '--', 'Color', colors{i-index+1}, 'Marker', 'x', 'DisplayName', 'Power Points')
    plot(regen_E_vec, eff_points_matrix(:,i,1) - maxi, '--', 'Color', colors{i-index+1}, 'Marker', 's', 'DisplayName', 'Efficiency Points')
    
end


xlabel('kWh Recovered')
ylabel('Points')
legend
grid on
improvePlot

figure;
hold all
for i = index:sweep_size
    plot(regen_E_vec, endurance_time_matrix(:,i,1), 'Color', colors{i-index+1}, 'DisplayName', ['Avg Accel Pwr = ' num2str(P_draw_no_regen_vec(i)/1000) ])
end

xlabel('kWh Recovered')
ylabel('Endurance Time')
legend
grid on
improvePlot

figure;
hold all
for i = index:sweep_size
    plot(regen_E_vec, t_to_60_matrix(:,i,1), 'Color', colors{i-index+1}, 'DisplayName', ['Avg Accel Pwr = ' num2str(P_draw_no_regen_vec(i)/1000) ])
end

xlabel('kWh Recovered')
ylabel('Time to 60C')
legend
grid on
improvePlot