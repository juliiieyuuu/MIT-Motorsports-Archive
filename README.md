# MIT-Motorsports-Archive
NOTE: The full battery sizing code is in a private repository that belongs to MIT Motorsports, so sharing access is restricted. For the sake of sharing a code sample, I have copied over a few scripts. Please refer to the summary below for a full explanation of the simulation and all the features it offers.

# MY21 Battery Sizing

This README will explain which scripts to run and describe how the SIM is organized. 

(See https://wikis.mit.edu/confluence/pages/viewpage.action?pageId=150209243 for a more conceptual discussion of BatteryPointsModel)

<b> How to Produce Plots: </b> \
The main script to run is BatteryPointsModel. Within this function, there are multiple sections. Depending on which plots you want to generate, the sections you run will vary. Here is a general breakdown:
- <b> Setup constants and input values + Choose which pack parameters you want: </b> This is where the user can change battery parameters (ex. Sensitivity, regen/no_regen favorable, lower/upper bound for pcount sweep, sweep size, battery pack, etc.) 
- <b> Best cell for each battery capacity: </b> Sweeps different battery capacities and outputs a graph showing the optimal cell for each capacity. Only run this section if you want to compare different cells. (See comments in script for details about calculating the reference point for cell comparison. I would recommend calculating the reference point once, then hard-coding it to save running time!)
- <b> Model Points: </b> This section calculates overall points, mass points, efficiency points, time to 60C, endurance time, number of unlimped laps, etc. using the parameters defined in the first section. (necessary for rainbow plots!) (See comments in script for details about calculating minimum laptime for a specific cell. This is needed to normalize all of our calculations in this section. I would recommend calculating this once, then hard-coding it to save running time!)
- <b> Make contour plots: </b> Interpolates the outputs from “Model Points” to generate rainbow contour plots! (Graphs each variable as a function of Average Power Draw)
- <b> Get best points for each pcount to compare pcounts: </b> Outputs a graph with pcount on the x axis and best points on the y axis (labels pcount that has the highest number of points)
- <b> Best Points vs P_draw_no_regen (for each p_count): </b> For each power draw value, the script outputs the highest points that can be attained by each pcount

<b> thermals_test_function: </b>
This script generates graphs that track power cool, power dissipation, distance traveled, battery temperature, battery voltage, etc. over time. Like BatteryPointsModel, certain battery parameters can be changed (located at the very top of the script)

<b> SIM Structure: </b>
- <b> Main: </b> BatteryPointsModel, thermals_test_function
- <b> Helper Functions: </b> find_point_parameters_new, run_endurance_thermals, thermals, find_laptime_from_power, split_power_or_current, delta_pts_battery, getAvgPw
