# call_fast_LP.jl
# Wrapper for fast_LP.jl, for Sloan computing clusters
# Authored by Arthur J Delarue on 10/28/15

include("intro.jl")
include("fast_LP.jl")

function load_parameters()
	if length(ARGS) == 0
		println("Error: shell script failed to pass in arguments")
		return nothing
	else
		JSONfileName = ARGS[1]
		params = JSON.parsefile(JSONfileName)
		return JSONfileName, params
	end
end

t=time()
# Load prameters
JSONfile, params = load_parameters()
# Load Manhattan
manhattan = loadCityGraph()
ws = params["WARMSTART"]
if ws
	startTimes = load("Outputs/fast_r200_minr1_i100_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st1000_add1000_every1_tc0.0/manhattan-times-100.jld", "times")
else
	startTimes = manhattan.roadTime
end
# Load data
if params["LOCALIZED_TESTING"]
	travel_times, num_rides, trainTestDf, testing_travel_times, testing_num_rides, testDf = loadInputTravelTimes(manhattan.positions, params["METHOD"], year=params["YEAR"], startTime=params["START_TIME"], endTime=params["END_TIME"], startMonth=params["START_MONTH"], endMonth=params["END_MONTH"], radius = params["RADIUS"], split_training_type = "localized")
else
	travel_times, num_rides, trainTestDf, testing_travel_times, testing_num_rides, testDf = loadInputTravelTimes(manhattan.positions, params["METHOD"], year=params["YEAR"], startTime=params["START_TIME"], endTime=params["END_TIME"], startMonth=params["START_MONTH"], endMonth=params["END_MONTH"], radius = params["RADIUS"], split_training_type = "random")
end
# Choose constraints if random
if params["RANDOM_CONSTRAINTS"]
	travel_times, num_rides = chooseConstraints(travel_times, num_rides, sample_size=params["SAMPLE_SIZE"]);
end
println("**** Launching computation ****")
# Launch ! (three cases for different error computations)
if params["ERROR_COMPUTATION"] == "single-ride"
	@time status, new_times, directory = fast_LP(manhattan, travel_times, num_rides, trainTestDf, startTimes, errorComputation=params["ERROR_COMPUTATION"], model_type=params["MODEL"], last_smooth=params["LAST_SMOOTH"], max_rounds=params["MAX_ROUNDS"], min_rides=params["MIN_RIDES"], radius=params["RADIUS"], year=params["YEAR"], startTime=params["START_TIME"], endTime=params["END_TIME"], startMonth=params["START_MONTH"], endMonth=params["END_MONTH"], method=params["METHOD"], sample_size=params["SAMPLE_SIZE"], turnCost=params["TURN_COST"], turnCostAsVariable=params["TURN_COST_AS_VARIABLE"], delta_bound=params["DELTA_BOUND"], maxNumPathsPerOD=params["MAX_NUM_PATHS_PER_OD"], computeFinalSP = params["COMPUTE_FINAL_SHORTEST_PATHS"], startWithSimpleLP=params["START_SIMPLE"], randomConstraints=params["RANDOM_CONSTRAINTS"], dynamicConstraints=params["DYNAMIC_CONSTRAINTS"], numPairsToAdd = params["NUM_OD_ADDED"], numPairsToRemove = params["NUM_OD_REMOVED"], addOnlyIfAbovePercentile=params["ADD_IF_ABOVE_PERCENTILE"], removeOnlyIfBelowPercentile=params["REMOVE_IF_BELOW_PERCENTILE"], globalConstraintUpdate=params["GLOBAL_CONSTRAINT_UPDATE"],iterationMultiple = params["UPDATE_EVERY_N_ITERATIONS"], warmstart=params["WARMSTART"], metropolis=params["METROPOLIS"], localized_testing=params["LOCALIZED_TESTING"])
elseif params["ERROR_COMPUTATION"] == "average"
	@time status, new_times, directory = fast_LP(manhattan, travel_times, num_rides, testing_travel_times, startTimes, errorComputation=params["ERROR_COMPUTATION"], model_type=params["MODEL"], last_smooth=params["LAST_SMOOTH"], max_rounds=params["MAX_ROUNDS"], min_rides=params["MIN_RIDES"], radius=params["RADIUS"], year=params["YEAR"], startTime=params["START_TIME"], endTime=params["END_TIME"], startMonth=params["START_MONTH"], endMonth=params["END_MONTH"], method=params["METHOD"], sample_size=params["SAMPLE_SIZE"], turnCost=params["TURN_COST"], turnCostAsVariable=params["TURN_COST_AS_VARIABLE"], delta_bound=params["DELTA_BOUND"], maxNumPathsPerOD=params["MAX_NUM_PATHS_PER_OD"], computeFinalSP = params["COMPUTE_FINAL_SHORTEST_PATHS"], startWithSimpleLP=params["START_SIMPLE"], randomConstraints=params["RANDOM_CONSTRAINTS"], dynamicConstraints=params["DYNAMIC_CONSTRAINTS"], numPairsToAdd = params["NUM_OD_ADDED"], numPairsToRemove = params["NUM_OD_REMOVED"], addOnlyIfAbovePercentile=params["ADD_IF_ABOVE_PERCENTILE"], removeOnlyIfBelowPercentile=params["REMOVE_IF_BELOW_PERCENTILE"], globalConstraintUpdate=params["GLOBAL_CONSTRAINT_UPDATE"],iterationMultiple = params["UPDATE_EVERY_N_ITERATIONS"], warmstart=params["WARMSTART"], metropolis=params["METROPOLIS"], localized_testing=params["LOCALIZED_TESTING"])
elseif params["ERROR_COMPUTATION"] == "both"
	@time status, new_times, directory = fast_LP(manhattan, travel_times, num_rides, trainTestDf, startTimes, testingData2 = testing_travel_times, errorComputation=params["ERROR_COMPUTATION"], model_type=params["MODEL"], last_smooth=params["LAST_SMOOTH"], max_rounds=params["MAX_ROUNDS"], min_rides=params["MIN_RIDES"], radius=params["RADIUS"], year=params["YEAR"], startTime=params["START_TIME"], endTime=params["END_TIME"], startMonth=params["START_MONTH"], endMonth=params["END_MONTH"], method=params["METHOD"], sample_size=params["SAMPLE_SIZE"], turnCost=params["TURN_COST"], turnCostAsVariable=params["TURN_COST_AS_VARIABLE"], delta_bound=params["DELTA_BOUND"], maxNumPathsPerOD=params["MAX_NUM_PATHS_PER_OD"], computeFinalSP = params["COMPUTE_FINAL_SHORTEST_PATHS"], startWithSimpleLP=params["START_SIMPLE"], randomConstraints=params["RANDOM_CONSTRAINTS"], dynamicConstraints=params["DYNAMIC_CONSTRAINTS"], numPairsToAdd = params["NUM_OD_ADDED"], numPairsToRemove = params["NUM_OD_REMOVED"], addOnlyIfAbovePercentile=params["ADD_IF_ABOVE_PERCENTILE"], removeOnlyIfBelowPercentile=params["REMOVE_IF_BELOW_PERCENTILE"], globalConstraintUpdate=params["GLOBAL_CONSTRAINT_UPDATE"],iterationMultiple = params["UPDATE_EVERY_N_ITERATIONS"], warmstart=params["WARMSTART"], metropolis=params["METROPOLIS"], localized_testing=params["LOCALIZED_TESTING"])
end
# Move log files
f=open(string(directory, "timestats.csv"), "a")
write(f, string(time()-t, '\n'))
close(f)
mv(JSONfile, string(directory, "parameters.json"), remove_destination=true)