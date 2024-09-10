classdef BatteryPointsClass < handle
    
    properties
        options_sim
        options_sweep
        battery
    end
    
    methods
        
        plot_points(obj)
                
        [points_index, points_matrix_best, mass_points_matrix_best, power_points_matrix_best, eff_points_matrix_best, ...
                t_to_60_matrix_best, p_count_matrix_best, endurance_time_matrix_best, unlimped_laps_matrix, P_draw_no_regen_vec] = ...
                sweep_run(obj)
                        
        [efficiency, eff_pts, Eyour] = eff_sens_run(obj)
            
        min_endur_t = find_min_endur_t_specific_cell(obj, battery, options_sweep, options_sim)
  
    end
    
    methods (Static)
      
        [points, mass_points, limp_power_points, eff_points] = delta_pts_battery(delta_mass_kg, endur_times, tmin, net_energies, Emin)
        
        [endur_times, E_used_total, P_draw_no_regen_vec, P_draw_limped_no_regen, time_to_60C, mass_delta, unlimped_laps, log] = find_point_parameters_transient(regen_E, eff, battery, options_sweep, options_sim)

         function obj = BatteryPointsClass(opt_struct)
           names = properties(obj);
           for i = 1:length(names)
               if ~isfield(opt_struct, names{i})
                   if isempty(obj.(names{i}))
                        error(strcat("'", names{i}, "' is not defined."))
                   end
               else
                    obj.(names{i}) = opt_struct.(names{i});
               end
           end
         end
        
    end    
    
end