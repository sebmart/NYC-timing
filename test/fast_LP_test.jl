# fast_LP_test.jl
# Interface for running fast_LP methods on metropolis instead of Manhattan
# This will call fast_LP with a modified version of Manhattan

# Default parameters: ignore since everything is loaded from JSON
PROB = 0.7
MODEL = "relaxed_smooth"
LAST_SMOOTH = false
MAX_ROUNDS = 50
MIN_RIDES = 1
MAX_NUM_PATHS_PER_OD = 3
TURN_COST = 2.0
TURN_COST_AS_VARIABLE = false
DELTA_BOUND = [0.6,0.5,0.4]
START_SIMPLE = false
DYNAMIC_CONSTRAINTS = true 			# true for dynamic constraints, false otw
GLOBAL_CONSTRAINT_UPDATE = false 	# true if we remove worst paths globally, false otw
SAMPLE_SIZE = 1000 					# starting number of constraints
NUM_OD_ADDED = 1000 				# number of (O,D) pairs to add
NUM_OD_REMOVED = 1000				# number of (O,D) pairs to remove
UPDATE_EVERY_N_ITERATIONS = 1 		# add number of (O,D) above every $N iterations
ADD_IF_ABOVE_PERCENTILE = 0.7		# only add rides that have an error above this percentile
REMOVE_IF_BELOW_PERCENTILE = 0.3	# only remove rides that have an error below this percentile

function fast_LP_test(;
	model_type::String=MODEL,
	last_smooth::Bool=LAST_SMOOTH,
	prob::Float64=PROB,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	turnCost::Float64=TURN_COST,
	turnCostAsVariable::Bool=TURN_COST_AS_VARIABLE,
	delta_bound::Union{Array{Float64,1},Array{Any,1}}=DELTA_BOUND,
	maxNumPathsPerOD::Int=MAX_NUM_PATHS_PER_OD,
	dynamicConstraints::Bool=DYNAMIC_CONSTRAINTS,
	numPairsToAdd::Int = NUM_OD_ADDED,
	numpairsToRemove::Int = NUM_OD_REMOVED,
	iterationMultiple::Int = UPDATE_EVERY_N_ITERATIONS,
	addOnlyIfAbovePercentile::Float64 = ADD_IF_ABOVE_PERCENTILE,
	removeOnlyIfBelowPercentile::Float64 = REMOVE_IF_BELOW_PERCENTILE,
	globalConstraintUpdate::Bool = GLOBAL_CONSTRAINT_UPDATE,
	sample_size::Int=SAMPLE_SIZE,
	startWithSimpleLP::Bool=START_SIMPLE)
	
	# This code loads the data
	adjacencyList = load("Inputs/input-graph.jld", "adjList")
	graph = DiGraph(length(adjacencyList))
	for i = 1:nv(graph), j in adjacencyList[i]
		add_edge!(graph, i, j)
	end
	real_times = load("Inputs/input-realTimes.jld", "realTimes")
	minRoadTime = load("Inputs/input-speedLimits.jld", "speedLimits")
	meanRoadTime = load("Inputs/input-meanTimes.jld", "meanTimes")
	distances = load("Inputs/input-distances.jld", "distances")
	coordinates = load("Inputs/input-positions.jld", "coordinates")
	positions = Coordinates[]
	for i = 1:nv(graph)
		push!(positions, Coordinates(coordinates[i,1], coordinates[i,2]))
	end

	travel_times = load("Inputs/input-travelTimes-$(prob).jld", "travelTimes")
	num_rides = load("Inputs/input-numRides-$(prob).jld", "numRides")

	metropolis_manhattan = Manhattan(emptyType = true)
	metropolis_manhattan.network = graph
	metropolis_manhattan.positions = positions
	metropolis_manhattan.distances = distances
	metropolis_manhattan.roadTime = minRoadTime

	status, new_times, directory = fast_LP(metropolis_manhattan, travel_times, num_rides, zeros(nv(graph),nv(graph)), minRoadTime, model_type=model_type, last_smooth = last_smooth, max_rounds=max_rounds, min_rides=min_rides, year=0, sample_size=sample_size, turnCost=turnCost, turnCostAsVariable=turnCostAsVariable, delta_bound=delta_bound, maxNumPathsPerOD=maxNumPathsPerOD, computeFinalSP=true,	startWithSimpleLP=startWithSimpleLP, randomConstraints=false, dynamicConstraints=dynamicConstraints, numPairsToAdd=numPairsToAdd, numPairsToRemove=numPairsToRemove, iterationMultiple = iterationMultiple,	addOnlyIfAbovePercentile = addOnlyIfAbovePercentile, removeOnlyIfBelowPercentile = removeOnlyIfBelowPercentile, globalConstraintUpdate=globalConstraintUpdate, warmstart=false, metropolis=true, real_TOD_metropolis = real_times, real_tij_metropolis = meanRoadTime, prob=prob, last_smooth=last_smooth, warmstart=false)

	return status, new_times, directory
end

# Load dependencies
include("metropolis_intro.jl")

t=time()
# Load parameters
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
JSONfile, params = load_parameters()

# run code
status, new_times, directory = @time fast_LP_test(model_type=params["MODEL"], last_smooth=params["LAST_SMOOTH"], prob=params["PROB"], max_rounds=params["MAX_ROUNDS"], min_rides=params["MIN_RIDES"], turnCost=params["TURN_COST"], turnCostAsVariable=params["TURN_COST_AS_VARIABLE"], delta_bound=params["DELTA_BOUND"], maxNumPathsPerOD=params["MAX_NUM_PATHS_PER_OD"], dynamicConstraints=params["DYNAMIC_CONSTRAINTS"], numPairsToAdd = params["NUM_OD_ADDED"], numpairsToRemove = params["NUM_OD_REMOVED"], iterationMultiple = params["UPDATE_EVERY_N_ITERATIONS"], addOnlyIfAbovePercentile = params["ADD_IF_ABOVE_PERCENTILE"], removeOnlyIfBelowPercentile = params["REMOVE_IF_BELOW_PERCENTILE"], globalConstraintUpdate = params["GLOBAL_CONSTRAINT_UPDATE"], sample_size=params["SAMPLE_SIZE"], startWithSimpleLP=params["START_SIMPLE"])

# Move log files
f=open(string(directory, "timestats.csv"), "a")
write(f, string(time()-t, '\n'))
close(f)
mv(JSONfile, string(directory, "parameters.json"), remove_destination=true)