# new_LP.jl
# New LP formulation for travel time estimation
# Authored by Arthur J Delarue on 7/2/15

MODEL = "relaxed_smooth"		# relaxed/strict _ smooth/nothing/random
LAST_SMOOTH = false 			# true if last iteration is smooth, false otherwise. Irrelevant if model is smooth already.
ERROR_COMPUTATION = "average" 		# average/single-ride/both

MAX_ROUNDS = 2					# max number of iterations

MIN_RIDES = 1					# min number of rides
RADIUS = 140					# radius (4D Euclidian)
YEAR = 2013						
START_TIME = 12
END_TIME = 14
START_MONTH = 1
END_MONTH = 1
WEEKDAYS = true 				# true for weekdays, false for weekends
METHOD = "average"				# average/nearest [neighbor]

RANDOM_CONSTRAINTS = false 		# true for purely random constraints, false otw
DYNAMIC_CONSTRAINTS = true 		# true for dynamic constraints, false otw
SAMPLE_SIZE = 1000 				# starting number of constraints
NUM_OD_ADDED = 1000 			# number of (O,D) pairs to add
UPDATE_EVERY_N_ITERATIONS = 1 	# add number of (O,D) above every $N iterations
SELECTION_RULE_CUTOFF = 0.9		# Value must be between 0. and 1., default is 0.9

TURN_COST = 0.	 				# turning cost initial value
TURN_COST_AS_VARIABLE = false 	# true if LP updates turning cost, false otw

DELTA_BOUND = [0.06, 0.05, 0.04]	# relaxation base (sketchy)

MAX_NUM_PATHS_PER_OD = 3 		# at least 1 please

COMPUTE_FINAL_SHORTEST_PATHS = true 	# Should always be true unless trying to save time somehow

START_SIMPLE = false 			# true if initial simple LP is used, false otw

METROPOLIS = false 				# always set this to false unless you wish to crash everything

function fast_LP(
	manhattan::Manhattan,
	travelTimes::Array{Float64,2},
	numRides::Array{Int,2},
	testingData::Union{DataFrame, Array{Float64,2}},
	startTimes::AbstractArray{Float64,2};
	errorComputation::AbstractString=ERROR_COMPUTATION,
	model_type::AbstractString=MODEL,
	last_smooth::Bool=LAST_SMOOTH,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	radius::Int=RADIUS,
	year::Int=YEAR,
	startTime::Int=START_TIME,
	endTime::Int=END_TIME,
	startMonth::Int=START_MONTH,
	endMonth::Int=END_MONTH,
	method::AbstractString=METHOD,
	sample_size::Int=SAMPLE_SIZE,
	turnCost::Float64=TURN_COST,
	turnCostAsVariable::Bool=TURN_COST_AS_VARIABLE,
	delta_bound::Union{Array{Float64,1},Array{Any,1}}=DELTA_BOUND,
	maxNumPathsPerOD::Int=MAX_NUM_PATHS_PER_OD,
	computeFinalSP::Bool = COMPUTE_FINAL_SHORTEST_PATHS,
	startWithSimpleLP::Bool=START_SIMPLE,
	randomConstraints::Bool=RANDOM_CONSTRAINTS,
	dynamicConstraints::Bool=DYNAMIC_CONSTRAINTS,
	numPairsToAdd::Int = NUM_OD_ADDED,
	iterationMultiple::Int = UPDATE_EVERY_N_ITERATIONS,
	selectionRuleCutoff::Float64 = SELECTION_RULE_CUTOFF,
	metropolis::Bool=METROPOLIS,
	real_TOD_metropolis::AbstractArray{Float64}=zeros(1,1),
	real_tij_metropolis::AbstractArray{Float64}=zeros(1,1),
	prob::Float64=0.0,
	testingData2::Array{Float64,2}=zeros(1,1))

	graph = manhattan.network
	roadTimes = deepcopy(manhattan.roadTime)
	distances = manhattan.distances
	positions = manhattan.positions
	roads = edges(graph)
	nodes = vertices(graph)
	out = [copy(out_neighbors(graph,i)) for i in nodes]
	inn = [copy(in_neighbors(graph,i)) for i in nodes]

	# Find nonzero travel times
	if DYNAMIC_CONSTRAINTS
		numDataPoints = sample_size
		pairs = [copy(Int[]) for i = 1:length(nodes)]
	else
		pairs = [find(travelTimes[i,:]) for i=nodes]
		numDataPoints = sum([length(pairs[i]) for i=nodes])
	end

	# Create output directory name (sorry this is so complicated)
	if !(last_smooth)
		TESTDIR = "fast_r$(radius)_minr$(min_rides)_i$(max_rounds)_wd_$(startTime)$(endTime)_m$(startMonth)$(endMonth)_y$(year)_$(model_type)_ppc$(maxNumPathsPerOD)_data$(method)_error$errorComputation"
	else
		TESTDIR = "fast_r$(radius)_minr$(min_rides)_i$(max_rounds)_wd_$(startTime)$(endTime)_m$(startMonth)$(endMonth)_y$(year)_$(model_type)_lsmooth_ppc$(maxNumPathsPerOD)_data$(method)_error$errorComputation"
	end
	if randomConstraints
		TESTDIR=string(TESTDIR, "_rnd_rides$(sample_size)")
	elseif dynamicConstraints
		TESTDIR=string(TESTDIR, "_dynConstr_st$(sample_size)_add$(numPairsToAdd)_every$(iterationMultiple)_sc$(selectionRuleCutoff)")
	end
	if turnCostAsVariable
		TESTDIR=string(TESTDIR, "_tcvar_start$(turnCost)")
	else
		TESTDIR=string(TESTDIR, "_tc$(turnCost)")
	end
	if startWithSimpleLP
		TESTDIR=string(TESTDIR, "_ss")
	end
	if metropolis
		TESTDIR=string(TESTDIR, "_metropolis$(prob)")
	end

	# Create directory if necessary:
	if !isdir("Outputs/$TESTDIR")
		mkdir("Outputs/$TESTDIR")
	end
	# Add vital information to file
	writeDataToFile("Outputs/$TESTDIR", model_type, max_rounds, min_rides, radius, false, 0, sample_size, turnCost)
	# Turn cost file
	if turnCostAsVariable
		turnCostFile = open("Outputs/$(TESTDIR)/tc.csv", "w")
	end
	# Metropolis error output
	if metropolis
		# specific error output
	else
		errorFile = open("Outputs/$(TESTDIR)/errorstats.csv", "w")
	end
	# Runtime info output
	timeFile = open("Outputs/$(TESTDIR)/timestats.csv", "w")
	if dynamicConstraints
		numConstraintsFile = open("Outputs/$(TESTDIR)/numODpairs.csv", "w")
	end

	# Tell us where output is going
	println("-- Saving outputs to directory Outputs/$(TESTDIR)/")

	# Create path storing array and initialize. By convention totalPaths[i,1] is the equality constraint
	numPaths = zeros(Int, numDataPoints)
	totalPaths = Array(Any, numDataPoints)
	totalNumExpensiveTurns = Array(Any, numDataPoints)
	for i = 1:numDataPoints
		totalPaths[i] = Array{Int}[]
		totalNumExpensiveTurns[i] = Int[]
		sizehint!(totalPaths[i], maxNumPathsPerOD)
		sizehint!(totalNumExpensiveTurns[i], maxNumPathsPerOD)
	end

	if errorComputation == "single-ride" || errorComputation == "both"
		println("**** Building KDTree of nodes for testing ****")
		nodeTree, nodePairs = buildNodeKDTree(manhattan.positions)
	end

	println("**** Initializing dataset ****")
	if dynamicConstraints
		srcs, dsts, newTravelTimes, newNumRides, pairs = chooseStartingConstraints(travelTimes, numRides, pairs, sample_size = sample_size, num_nodes = length(nodes))
	else
		# Run over all pairs of nodes that have data
		srcs = Int[]
		dsts = Int[]
		sizehint!(srcs, sample_size)
		sizehint!(dsts, sample_size)
		for i in nodes, j in nodes
			if travelTimes[i,j] > 0
				# Load sources and destinations: will be useful when adding constraints
				push!(srcs, i)
				push!(dsts, j)
			end
		end
	end

	status = 0
	newTimes = startTimes
	old_objective = 0

	# Compute shortest paths (with turn cost)
	if startWithSimpleLP
		println("**** Finding initial set of paths ****")
		status, newTimes, turnCost, totalPaths, totalNumExpensiveTurns = simple_LP(manhattan, travelTimes, numRides, testingData, TESTDIR, model_type=split(model_type, "_")[1], max_rounds=max_rounds, turnCost=turnCost, turnCostAsVariable=turnCostAsVariable, delta_bound=delta_bound, maxNumPathsPerOD=maxNumPathsPerOD, computeFinalSP=false)
	else
		println("**** Computing shortest paths ****")
		new_graph, new_edge_dists, new_nodes, new_edge_isExpensive = modifyGraphForDijkstra(graph, newTimes, positions, turn_cost=turnCost)
		old_nodes = getInverseMapping(new_nodes, nv(new_graph))
		@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
	end

	# l is round number
	l = 1
	while l <= max_rounds
		println("###### ROUND $l ######")
		tic()
		actual_l = l
		# Update bound on relaxation variables
		if split(model_type, "_")[1] == "relaxed"
			delta = (1/2)^(div(actual_l-1,length(delta_bound)))*delta_bound[((actual_l-1)%length(delta_bound)) + 1]
			# Avoid numerical issues
			if delta < 1e-4 || l == max_rounds
				delta = 0
			end
		end

		# Create JuMP model for LP and QP
		println("**** Creating LP instance ****")
		if l == 1 && !startWithSimpleLP && !dynamicConstraints
			tolerance = 1e-4
		else
			tolerance = 1e-6
		end
		m = Model(solver=GurobiSolver(TimeLimit=10000, OutputFlag=1, Method=3, BarConvTol=tolerance))

		# Add one variable for each road, for each model
		@defVar(m, t[i=nodes,j=out[i]] >= 0.038 * manhattan.distances[i,j])

		# Add decision variable for turn cost if necessary
		if turnCostAsVariable
			@defVar(m, tc >= 0)
		else
			tc = turnCost
		end

		# Epsilon is the variable used for the absolute value difference between T and \hat{T}
		@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)

		# Create handles for inequality constraints
		if maxNumPathsPerOD > 1
			@defConstrRef inequalityPath[1:numDataPoints, 1:(max_rounds-1)]
		end

		# Set first LP objective
		@setObjective(m, Min, sum{ sqrt(numRides[i,j]/travelTimes[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]})

		# Add path constraints
		println("**** Adding constraints ****")
		if !(startWithSimpleLP) || l > 1
			println("---Fetching paths")
			paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurns(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, new_edge_isExpensive)
			println("---Updating path array")
		end
		for i=1:numDataPoints
			if !(startWithSimpleLP) || l > 1
				index = findfirst(totalPaths[i], paths[i])
				# If path already in paths
				if index != 0
					totalPaths[i][1], totalPaths[i][index] = totalPaths[i][index], totalPaths[i][1]
					totalNumExpensiveTurns[i][1], totalNumExpensiveTurns[i][index] = totalNumExpensiveTurns[i][index], totalNumExpensiveTurns[i][1]			# New path is equality constraint
				# If we still have space to add the path
				elseif numPaths[i] < maxNumPathsPerOD
					numPaths[i] += 1
					if numPaths[i] == 1
						push!(totalPaths[i], paths[i])
						push!(totalNumExpensiveTurns[i], numExpensiveTurns[i])
					else
						push!(totalPaths[i], paths[i])
						push!(totalNumExpensiveTurns[i], numExpensiveTurns[i])
						assert(numPaths[i] == length(totalPaths[i]))
						totalPaths[i][numPaths[i]], totalPaths[i][1] = totalPaths[i][1], totalPaths[i][numPaths[i]]
						totalNumExpensiveTurns[i][numPaths[i]], totalNumExpensiveTurns[i][1] = totalNumExpensiveTurns[i][1], totalNumExpensiveTurns[i][numPaths[i]]
					end
				# If we need to remove a path
				else
					worstIndex = findWorstPathIndex(totalPaths[i], totalNumExpensiveTurns[i], turnCost, newTimes)
					if worstIndex == 1
						totalPaths[i][1] = paths[i]
						totalNumExpensiveTurns[i][1] = numExpensiveTurns[i]
					else
						totalPaths[i][1], totalPaths[i][worstIndex] = paths[i], totalPaths[i][1]
						totalNumExpensiveTurns[i][1], totalNumExpensiveTurns[i][worstIndex] = numExpensiveTurns[i], totalNumExpensiveTurns[i][1]
					end
				end
			else
				numPaths[i] = length(totalPaths[i])
			end
			# Add inequality constraints
			if numPaths[i] > 1
				for j = 1:(numPaths[i]-1)
					if split(model_type, "_")[1] == "strict"
						inequalityPath[i,j] = @addConstraint(m, sum{t[totalPaths[i][j+1][a],totalPaths[i][j+1][a+1]], a=1:(length(totalPaths[i][j+1])-1)} - sum{t[totalPaths[i][1][a],totalPaths[i][1][a+1]], a=1:(length(totalPaths[i][1])-1)} >= - tc * totalNumExpensiveTurns[i][j+1] + tc * totalNumExpensiveTurns[i][1])
					elseif split(model_type, "_")[1] == "relaxed"
						inequalityPath[i,j] = @addConstraint(m, sum{t[totalPaths[i][j+1][a],totalPaths[i][j+1][a+1]], a=1:(length(totalPaths[i][j+1])-1)} - sum{t[totalPaths[i][1][a],totalPaths[i][1][a+1]], a=1:(length(totalPaths[i][1])-1)} >= - tc * totalNumExpensiveTurns[i][j+1] + tc * totalNumExpensiveTurns[i][1] - delta * travelTimes[srcs[i],dsts[i]])
					end
				end
			end
		end
		# Equality constraints (shortest path close to travel time data)
		@addConstraint(m, TLowerBound[i=1:numDataPoints], sum{t[totalPaths[i][1][a],totalPaths[i][1][a+1]], a=1:(length(totalPaths[i][1])-1)} + tc * totalNumExpensiveTurns[i][1] - travelTimes[srcs[i],dsts[i]] >= - epsilon[srcs[i],dsts[i]])
		@addConstraint(m, TUpperBound[i=1:numDataPoints], sum{t[totalPaths[i][1][a],totalPaths[i][1][a+1]], a=1:(length(totalPaths[i][1])-1)} + tc * totalNumExpensiveTurns[i][1] - travelTimes[srcs[i],dsts[i]] <= epsilon[srcs[i],dsts[i]])

		# Define objective variables and constraints for second part of m
		if split(model_type, "_")[2] == "smooth" || (last_smooth && l == max_rounds)
			@defVar(m, delta2[i=nodes,j=out[i]] >= 0)
			@addConstraint(m, objConstrLower[i=nodes,j=out[i]], -1 * t[i,j]/distances[i,j] + 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t[j,k], k = out[j]} + sum{1/distances[h,i] * t[h,i], h=inn[i]}) <= delta2[i,j])
			@addConstraint(m, objConstrUpper[i=nodes,j=out[i]], t[i,j]/distances[i,j] - 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t[j,k], k = out[j]} + sum{1/distances[h,i] * t[h,i], h=inn[i]}) <= delta2[i,j])
		elseif split(model_type, "_")[2] == "random"
			coeffVec = 2 * rand(length(roads)) - 1
			coeffs = spzeros(length(nodes), length(nodes))
			k = 1
			for i = nodes, j = out[i]
				coeffs[i,j] = coeffVec[k]
				k += 1
			end
			if turnCostAsVariable
				coeff_tc = 2 * rand() - 1
			end
		end
		# Solve LP
		println("**** Solving LP ****")
		status = solve(m)
		if status == :Infeasible
			buildInternalModel(m)
			print_iis_gurobi(m)
			break
		end
		objective = getObjectiveValue(m)
		if turnCostAsVariable
			println("Left turn cost: ", getValue(tc))
		end

		println("**** Setting up second LP ****")
		if split(model_type, "_")[2] == "smooth" || (last_smooth && l == max_rounds)
			if turnCostAsVariable
				@setObjective(m, Min, sum{delta2[i,j], i=nodes, j=out[i]} + tc)
			else
				@setObjective(m, Min, sum{delta2[i,j], i=nodes, j=out[i]})
			end
		elseif split(model_type, "_")[2] == "random"
			# Need to minimize absolute value (otherwise unbounded) so define new variable
			@defVar(m, zAbsVal >= 0)
			# Depending on whether turn cost is a variable, add it (or not) to the objective
			if turnCostAsVariable
				@addConstraint(m, randomObjectiveConstraint1, sum{coeffs[i,j] * t[i,j]/distances[i,j], i=nodes, j=out[i]} + coeff_tc * tc <= zAbsVal)
				@addConstraint(m, randomObjectiveConstraint2, sum{coeffs[i,j] * t[i,j]/distances[i,j], i=nodes, j=out[i]} + coeff_tc * tc >= - zAbsVal)
			else
				@addConstraint(m, randomObjectiveConstraint1, sum{coeffs[i,j] * t[i,j]/distances[i,j], i=nodes, j=out[i]} <= zAbsVal)
				@addConstraint(m, randomObjectiveConstraint2, sum{coeffs[i,j] * t[i,j]/distances[i,j], i=nodes, j=out[i]} >= - zAbsVal)
			end		
			@setObjective(m, Min, zAbsVal)		
		end

		if split(model_type, "_")[2] == "nothing" && !(last_smooth && l == max_rounds)
			# do nothing
			if turnCostAsVariable
				turnCost = getValue(tc)
				write(turnCostFile, string(l,",", turnCost, "\n"))
			end
		else
			# Solve second LP
			println("**** Solving second LP ****")
			@addConstraint(m, fixObjective, sum{sqrt(numRides[i,j]/travelTimes[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]} <=objective + 1e-1)
			status = solve(m)
			# Get new turning cost
			if turnCostAsVariable
				turnCost = getValue(tc)
				println(turnCost)
				write(turnCostFile, string(l,",", turnCost, "\n"))
			end
		end

		# Debug if infeasible
		if status == :Infeasible
			println("!!!! Computing IIS !!!!")
			print_iis_gurobi(m)
			break
		# Prepare output
		elseif status == :Optimal || status == :Suboptimal
			result = getValue(t)
			newTimes = spzeros(length(nodes), length(nodes))
			for element in result
				newTimes[element[1], element[2]] = element[3]
			end
			# Save updated Manhattan road times to file
			if metropolis
				saveRoadTimes(newTimes, "$TESTDIR/metropolis-times-$l")
			else
				saveRoadTimes(newTimes, "$TESTDIR/manhattan-times-$l")
			end
			if abs(objective - old_objective)/old_objective < 1e-10
				save("Outputs/$TESTDIR/end.jld", "num_iter", l)
				l = max_rounds
			else
				old_objective = objective
			end
		elseif status == :UserLimit
			println("!!!! User time limit exceeded !!!!")
			break
		end
		if l < max_rounds || computeFinalSP
			println("**** Computing shortest paths ****")
			new_graph, new_edge_dists, new_nodes, new_edge_isExpensive = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
			old_nodes = getInverseMapping(new_nodes, nv(new_graph))
			@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
			if metropolis
				printErrorStatsToFile("Outputs/$TESTDIR/errorstats-$(actual_l).txt", real_TOD_metropolis, travelTimes, new_sp.traveltime, numRides, real_tij_metropolis, newTimes, num_nodes = nv(graph))
			else
				if errorComputation == "single-ride"
					@time avg_sq_err, avg_rel_err, avg_bias = computeTestingError(testingData, nodeTree, nodePairs, new_sp.traveltime, method)
					write(errorFile, string(l, ",", avg_sq_err, ",", avg_rel_err, ",", avg_bias,"\n"))
					flush(errorFile)
				elseif errorComputation == "average"
					@time avg_sq_err, avg_rel_err, avg_bias = computeAverageTestingError(testingData, new_sp.traveltime, num_nodes=nv(graph))
					write(errorFile, string(l, ",", avg_sq_err, ",", avg_rel_err, ",", avg_bias,"\n"))
					flush(errorFile)
				elseif errorComputation == "both"
					@time avg_sq_err, avg_rel_err, avg_bias = computeTestingError(testingData, nodeTree, nodePairs, new_sp.traveltime, method)
					@time avg_sq_err_2, avg_rel_err_2, avg_bias_2 = computeAverageTestingError(testingData2, new_sp.traveltime, num_nodes=nv(graph))
					write(errorFile, string(l, ",", avg_sq_err, ",", avg_rel_err, ",", avg_bias, ",", avg_sq_err_2, ",", avg_rel_err_2, ",", avg_bias_2,"\n"))
					flush(errorFile)
				end
			end
			if dynamicConstraints
				write(numConstraintsFile, string(l,",",length(srcs), "\n"))
				flush(numConstraintsFile)
				if l < max_rounds && l % iterationMultiple == 0
					println("**** Updating constraint set ****")
					srcs, dsts, totalPaths, totalNumExpensiveTurns, numPaths, pairs = updateConstraints(travelTimes, numRides, new_sp.traveltime, totalPaths, totalNumExpensiveTurns, numPaths, srcs, dsts, pairs, numNodePairsToAdd = numPairsToAdd, selectionRuleCutoff=selectionRuleCutoff)
					numDataPoints = length(srcs)
				end
			end
		end
		write(timeFile, string(l, ",", toc(), "\n"))
		flush(timeFile)
		l += 1
	end
	println("-- Saved outputs to directory Outputs/$(TESTDIR)/")
	# Close log files
	if turnCostAsVariable
		close(turnCostFile)
	end
	if !metropolis
		close(errorFile)
	end
	close(timeFile)
	if dynamicConstraints
		close(numConstraintsFile)
	end
	return status, newTimes, "Outputs/$TESTDIR/"
end

# manhattan = loadCityGraph()
# travel_times, num_rides, trainTestDf, testing_travel_times, testing_num_rides, testDf = loadInputTravelTimes(manhattan.positions, METHOD, year=YEAR, startTime=START_TIME, endTime=END_TIME, startMonth=START_MONTH, endMonth=END_MONTH, radius = RADIUS)
# if RANDOM_CONSTRAINTS
# 	travel_times, num_rides = chooseConstraints(travel_times, num_rides, sample_size=SAMPLE_SIZE);
# end
# if ERROR_COMPUTATION == "single-ride"
# 	@time status, new_times = fast_LP(manhattan, travel_times, num_rides, trainTestDf, manhattan.roadTime)
# elseif ERROR_COMPUTATION == "average"
# 	@time status, new_times = fast_LP(manhattan, travel_times, num_rides, testing_travel_times, manhattan.roadTime)
# elseif ERROR_COMPUTATION == "both"
# 	@time status, new_times = fast_LP(manhattan, travel_times, num_rides, trainTestDf, manhattan.roadTime, testingData2 = testing_travel_times)
# end