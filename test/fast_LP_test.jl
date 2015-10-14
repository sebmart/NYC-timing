# fast_LP_test.jl
# Interface for running fast_LP methods on metropolis instead of Manhattan
# This will call fast_LP with a modified version of Manhattan

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

DYNAMIC_CONSTRAINTS = true
NUM_OD_ADDED = 1000
SAMPLE_SIZE = 1500
UPDATE_EVERY_N_ITERATIONS = 5

function fast_LP_test(;
	model_type::String=MODEL,
	last_smooth::Bool=LAST_SMOOTH,
	prob::Float64=PROB,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	turnCost::Float64=TURN_COST,
	turnCostAsVariable::Bool=TURN_COST_AS_VARIABLE,
	delta_bound::Array{Float64}=DELTA_BOUND,
	maxNumPathsPerOD::Int=MAX_NUM_PATHS_PER_OD,
	dynamicConstraints::Bool=DYNAMIC_CONSTRAINTS,
	numPairsToAdd::Int = NUM_OD_ADDED,
	iterationMultiple::Int = UPDATE_EVERY_N_ITERATIONS,
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

	status, new_times = @time fast_LP(metropolis_manhattan, travel_times, num_rides, zeros(nv(graph),nv(graph)), minRoadTime, model_type=model_type, last_smooth = false, max_rounds = max_rounds, min_rides = min_rides, times = "0000", preprocess = false, turnCost = turnCost, turnCostAsVariable = turnCostAsVariable, delta_bound = delta_bound, maxNumPathsPerOD = maxNumPathsPerOD, computeFinalSP = true, startWithSimpleLP = startWithSimpleLP, randomConstraints=false, dynamicConstraints = dynamicConstraints, sample_size=sample_size, numPairsToAdd=numPairsToAdd, iterationMultiple=iterationMultiple, metropolis = true, real_TOD_metropolis = real_times, real_tij_metropolis = meanRoadTime, prob=prob, last_smooth=last_smooth)

	return status, new_times
end

status, new_times = fast_LP_test()

# turn_costs = [0.0, 1.0, 2.0, 3.0, 4.0, 10.0, 20.0]
# turn_costs_as_variable = [true, false]
# max_ppc = [3, 5, 10]
# for tc in turn_costs, tcvar in turn_costs_as_variable, mppc in max_ppc
# 	status, new_times = fast_LP_test(turnCost = tc, turnCostAsVariable = tcvar, maxNumPathsPerOD = mppc)
# end
