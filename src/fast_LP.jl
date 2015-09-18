# new_LP.jl
# New LP formulation for travel time estimation
# Authored by Arthur J Delarue on 7/2/15

MODEL = "relaxed_smooth"
MAX_ROUNDS = 10
MIN_RIDES = 1
RADIUS = 140
TIMES = "1214"

PREPROCESS = true
NUM_CLUSTERS = 50

SAMPLE_SIZE = 50000

RANDOM_CONSTRAINTS = false

TURN_COST = 20.
TURN_COST_AS_VARIABLE = false

DELTA_BOUND = [0.06, 0.05, 0.04]

MAX_NUM_PATHS_PER_OD = 3 # at least 1 please

COMPUTE_FINAL_SHORTEST_PATHS = true

START_SIMPLE = true

function fast_LP(
	manhattan::Manhattan,
	travelTimes::Array{Float64,2},
	numRides::Array{Int,2},
	testingData::Array{Float64,2},
	startTimes::AbstractArray{Float64,2};
	model_type::String=MODEL,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	radius::Int=RADIUS,
	times::String=TIMES,
	preprocess::Bool=PREPROCESS,
	num_clusters::Int=NUM_CLUSTERS,
	sample_size::Int=SAMPLE_SIZE,
	turnCost::Float64=TURN_COST,
	turnCostAsVariable::Bool=TURN_COST_AS_VARIABLE,
	delta_bound::Array{Float64}=DELTA_BOUND,
	maxNumPathsPerOD::Int=MAX_NUM_PATHS_PER_OD,
	computeFinalSP::Bool = COMPUTE_FINAL_SHORTEST_PATHS,
	startWithSimpleLP::Bool=START_SIMPLE,
	randomConstraints::Bool=RANDOM_CONSTRAINTS)

	graph = manhattan.network
	roadTimes = deepcopy(manhattan.roadTime)
	distances = manhattan.distances
	positions = manhattan.positions
	roads = edges(graph)
	nodes = vertices(graph)
	out = [copy(out_neighbors(graph,i)) for i in nodes]
	inn = [copy(in_neighbors(graph,i)) for i in nodes]

	# Find nonzero travel times
	pairs = [find(travel_times[i,:]) for i=nodes]
	numDataPoints = sum([length(pairs[i]) for i=nodes])

	TESTDIR = "fast_r$(radius)_minr$(min_rides)_i$(max_rounds)_wd_$(times)_$(model_type)_ppc$(maxNumPathsPerOD)"
	if preprocess
		TESTDIR=string(TESTDIR, "_clust$(num_clusters)_rides$(sample_size)")
	elseif randomConstraints
		TESTDIR=string(TESTDIR, "_rnd_rides$(sample_size)")
	end
	if turnCostAsVariable
		TESTDIR=string(TESTDIR, "_tcvar")
	else
		TESTDIR=string(TESTDIR, "_tc$(turnCost)")
	end
	if startWithSimpleLP
		TESTDIR=string(TESTDIR, "_ss")
	end

	# Create directory if necessary:
	if !isdir("Outputs/$TESTDIR")
		mkdir("Outputs/$TESTDIR")
	end
	# Add vital information to file
	writeDataToFile("Outputs/$TESTDIR", model_type, max_rounds, min_rides, radius, preprocess, num_clusters, sample_size, turnCost)
	if turnCostAsVariable
		turnCostFile = open("Outputs/$(TESTDIR)/tc.csv", "w")
	end
	errorFile = open("Outputs/$(TESTDIR)/errorstats.csv", "w")
	timeFile = open("Outputs/$(TESTDIR)/timestats.csv", "w")

	# Tell us where output is going
	println("-- Saving outputs to directory Outputs/$(TESTDIR)/")

	# Create path storing array and initialize. By convention totalPaths[i,1] is the equality constraint
	numPaths = zeros(Int, numDataPoints)
	totalPaths = Array(Array{Int},(numDataPoints, maxNumPathsPerOD))
	for i = 1:numDataPoints, j = 1:maxNumPathsPerOD
		totalPaths[i,j] = Int[]
	end
	totalNumExpensiveTurns = Array(Int, (numDataPoints, maxNumPathsPerOD))

	# Run over all pairs of nodes that have data
	srcs = Int[]
	dsts = Int[]
	sizehint(srcs, sample_size)
	sizehint(dsts, sample_size)
	for i in nodes, j in nodes
		if travelTimes[i,j] > 0
			# Load sources and destinations: will be useful when adding constraints
			push!(srcs, i)
			push!(dsts, j)
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
		@time new_graph, new_edge_dists, new_nodes, new_edge_isExpensive = modifyGraphForDijkstra(graph, newTimes, positions, turn_cost=turnCost)
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
		if l == 1 && !startWithSimpleLP
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
			@defConstrRef inequalityPath[1:numDataPoints, 1:(maxNumPathsPerOD-1)]
		end

		# Set first LP objective
		@setObjective(m, Min, sum{ sqrt(numRides[i,j]/travelTimes[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]})

		# Add path constraints
		println("**** Adding constraints ****")
		if !(startWithSimpleLP) || l > 1
			paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, new_edge_isExpensive)
		end
		for i=1:numDataPoints
			if !(startWithSimpleLP) || l > 1
				index = findfirst(totalPaths[i, 1:maxNumPathsPerOD], paths[i])
				# If path already in paths
				if index != 0
					totalPaths[i,1], totalPaths[i,index] = totalPaths[i,index], totalPaths[i,1]
					totalNumExpensiveTurns[i,1], totalNumExpensiveTurns[i,index] = totalNumExpensiveTurns[i,index], totalNumExpensiveTurns[i,1]			# New path is equality constraint
				# If we still have space to add the path
				elseif numPaths[i] < maxNumPathsPerOD
					numPaths[i] += 1
					if numPaths[i] == 1
						totalPaths[i,1] = paths[i]
						totalNumExpensiveTurns[i,1] = numExpensiveTurns[i]
					else
						totalPaths[i,numPaths[i]], totalPaths[i,1] = totalPaths[i,1], paths[i]
						totalNumExpensiveTurns[i,numPaths[i]], totalNumExpensiveTurns[i,1] = totalNumExpensiveTurns[i,1], numExpensiveTurns[i]
					end
				# If we need to remove a path
				else
					worstIndex = findWorstPathIndex(totalPaths[i,1:maxNumPathsPerOD], totalNumExpensiveTurns[i,1:maxNumPathsPerOD], turnCost, newTimes)
					if worstIndex == 1
						totalPaths[i,1] = paths[i]
						totalNumExpensiveTurns[i,1] = numExpensiveTurns[i]
					else
						totalPaths[i,1], totalPaths[i, worstIndex] = paths[i], totalPaths[i,1]
						totalNumExpensiveTurns[i,1], totalNumExpensiveTurns[i,worstIndex] = numExpensiveTurns[i], totalNumExpensiveTurns[i,1]
					end
				end
			else
				j = 1
				while j < maxNumPathsPerOD
					if length(totalPaths[i,j+1]) > 0
						j += 1
					else
						break
					end
				end
				numPaths[i] = j
			end
			# Add inequality constraints
			if numPaths[i] > 1
				for j = 1:(numPaths[i]-1)
					if split(model_type, "_")[1] == "strict"
						inequalityPath[i,j] = @addConstraint(m, sum{t[totalPaths[i,j+1][a],totalPaths[i,j+1][a+1]], a=1:(length(totalPaths[i,j+1])-1)} - sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} >= - tc * totalNumExpensiveTurns[i,j+1] + tc * totalNumExpensiveTurns[i,1])
					elseif split(model_type, "_")[1] == "relaxed"
						inequalityPath[i,j] = @addConstraint(m, sum{t[totalPaths[i,j+1][a],totalPaths[i,j+1][a+1]], a=1:(length(totalPaths[i,j+1])-1)} - sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} >= - tc * totalNumExpensiveTurns[i,j+1] + tc * totalNumExpensiveTurns[i,1] - delta * travelTimes[srcs[i],dsts[i]])
					end
				end
			end
		end
		# Equality constraints (shortest path close to travel time data)
		@addConstraint(m, TLowerBound[i=1:numDataPoints], sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} + tc * totalNumExpensiveTurns[i,1] - travelTimes[srcs[i],dsts[i]] >= - epsilon[srcs[i],dsts[i]])
		@addConstraint(m, TUpperBound[i=1:numDataPoints], sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} + tc * totalNumExpensiveTurns[i,1] - travelTimes[srcs[i],dsts[i]] <= epsilon[srcs[i],dsts[i]])

		# Define objective variables and constraints for second part of m
		if split(model_type, "_")[2] == "smooth"
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
			println(getValue(tc))
		end

		println("**** Setting up second LP ****")
		if split(model_type, "_")[2] == "smooth"
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

		if split(model_type, "_")[2] == "nothing"
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
			saveRoadTimes(newTimes, "$TESTDIR/manhattan-times-$l")
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
			@time new_graph, new_edge_dists, new_nodes, new_edge_isExpensive = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
			old_nodes = getInverseMapping(new_nodes, nv(new_graph))
			@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
			avg_sq_err, avg_rel_err, avg_bias = compute_error(testingData, new_sp.traveltime)
			write(errorFile, string(l, ",", avg_sq_err, ",", avg_rel_err, ",", avg_bias,"\n"))
		end
		l += 1
		write(timeFile, string(l, ",", toc(), "\n"))
	end
	println("-- Saved outputs to directory Outputs/$(TESTDIR)/")
	if turnCostAsVariable
		close(turnCostFile)
	end
	close(errorFile)
	close(timeFile)
	return status, newTimes
end

manhattan = loadCityGraph()
if PREPROCESS
	@time outputPreprocessedConstraints(manhattan, "training", radius=RADIUS, numClusters=NUM_CLUSTERS, minRides=MIN_RIDES, sampleSize=SAMPLE_SIZE, overwrite = false, times=TIMES)
end
if RANDOM_CONSTRAINTS
	travel_times, num_rides = loadNewTravelTimeData(trainOrTest = "training", radius=RADIUS, times=TIMES, min_rides=MIN_RIDES, preprocess=false, loadTestingMatrixDirectly=true)
	travel_times, num_rides = chooseConstraints(travel_times, num_rides, sample_size=SAMPLE_SIZE);
else
	travel_times, num_rides = loadNewTravelTimeData(trainOrTest = "training", radius=RADIUS, times=TIMES, min_rides=MIN_RIDES, preprocess=PREPROCESS, num_clusters=NUM_CLUSTERS, sampleSize = SAMPLE_SIZE);
end
testing_data, numRides = loadNewTravelTimeData(trainOrTest="training", radius = RADIUS, times = TIMES, preprocess = false, loadTestingMatrixDirectly = true, saveTestingMatrix = false);
@time status, new_times = fast_LP(manhattan, travel_times, num_rides, testing_data, manhattan.roadTime)