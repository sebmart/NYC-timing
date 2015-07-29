# new_LP.jl
# New LP formulation for travel time estimation
# Authored by Arthur J Delarue on 7/2/15

MODEL = "relaxed"
MAX_ROUNDS = 30
MIN_RIDES = 1
RADIUS = 140
TIMES = "1214"

PREPROCESS = true
NUM_CLUSTERS = 50
SAMPLE_SIZE = 5000

TURN_COST = 10.

DELTA_BOUND = [0.06,0.05,0.04]

MAX_NUM_PATHS_PER_OD = 3 # at least 1 please

function fast_LP(
	manhattan::Manhattan,
	travelTimes::Array{Float64,2},
	numRides::Array{Int,2},
	testingData::Array{Float64,2};
	model_type::String=MODEL,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	radius::Int=RADIUS,
	times::String=TIMES,
	preprocess::Bool=PREPROCESS,
	num_clusters::Int=NUM_CLUSTERS,
	sample_size::Int=SAMPLE_SIZE,
	turnCost::Float64=TURN_COST,
	delta_bound::Array{Float64}=DELTA_BOUND,
	maxNumPathsPerOD::Int=MAX_NUM_PATHS_PER_OD)

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

	TMAX = 3600
	TESTDIR = "fast_r$(radius)_minr$(min_rides)_i$(max_rounds)_wd_$(times)_$(model_type)_ppc$(maxNumPathsPerOD)"
	if preprocess
		TESTDIR=string(TESTDIR, "_clust$(num_clusters)_rides$(sample_size)")
	end

	# Create directory if necessary:
	if !isdir("Outputs/$TESTDIR")
		mkdir("Outputs/$TESTDIR")
	end
	# Add vital information to file
	writeDataToFile("Outputs/$TESTDIR", model_type, max_rounds, min_rides, radius, preprocess, num_clusters, sample_size, turnCost)

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
	newTimes = roadTimes
	old_objective = 0

	# Compute shortest paths (with turn cost)
	println("**** Computing shortest paths ****")
	@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, positions, turn_cost=turnCost)
	old_nodes = getInverseMapping(new_nodes, nv(new_graph))
	@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)

	# l is round number
	l = 1
	while l <= max_rounds
		println("###### ROUND $l ######")
		actual_l = l
		# Update bound on relaxation variables
		if model_type == "relaxed"
			delta = (1/2)^(div(actual_l-1,length(delta_bound)))*delta_bound[((actual_l-1)%length(delta_bound)) + 1]
			# Avoid numerical issues
			if delta < 1e-4 || l == max_rounds
				delta = 0
			end
		end

		# Create JuMP model for LP and QP
		println("**** Creating LP instance ****")
		m = Model(solver=GurobiSolver(TimeLimit=10000, OutputFlag=1, Method=3))

		# Add one variable for each road, for each model
		@defVar(m, t[i=nodes,j=out[i]] >= roadTimes[i,j])

		# Epsilon is the variable used for the absolute value difference between T and \hat{T}
		@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)
		# T is \hat{T} in the notes
		@defVar(m, T[i=nodes, j=pairs[i]] >= 0)
		@addConstraint(m, TLowerBound[i=nodes,j=pairs[i]], T[i,j] - travelTimes[i,j] >= - epsilon[i,j])
		@addConstraint(m, TUpperBound[i=nodes,j=pairs[i]], T[i,j] - travelTimes[i,j] <= epsilon[i,j])

		# Define objective variables and constraints for m (second part)
		@defVar(m, delta2[i=nodes,j=out[i]] >= 0)
		@addConstraint(m, objConstrLower[i=nodes,j=out[i]], -1 * t[i,j]/distances[i,j] + 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t[j,k], k = out[j]} + sum{1/distances[h,i] * t[h,i], h=inn[i]}) <= delta2[i,j])
		@addConstraint(m, objConstrUpper[i=nodes,j=out[i]], t[i,j]/distances[i,j] - 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t[j,k], k = out[j]} + sum{1/distances[h,i] * t[h,i], h=inn[i]}) <= delta2[i,j])

		# Create handles for inequality constraints
		if maxNumPathsPerOD > 1
			@defConstrRef inequalityPath[1:numDataPoints, 1:(maxNumPathsPerOD-1)]
		end

		# Set first LP objective
		@setObjective(m, Min, sum{ sqrt(numRides[i,j]/travelTimes[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]})

		# Add path constraints
		println("**** Adding constraints ****")
		tic()
		paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, newTimes, new_edge_dists)
		for i=1:numDataPoints
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
			# Add inequality constraints
			if numPaths[i] > 1
				for j = 1:(numPaths[i]-1)
					if model_type == "strict"
						inequalityPath[i,j] = @addConstraint(m, sum{t[totalPaths[i,j][a],totalPaths[i,j][a+1]], a=1:(length(totalPaths[i,j])-1)} - T[srcs[i],dsts[i]] >= - turnCost * totalNumExpensiveTurns[i,j])
					elseif model_type == "relaxed"
						inequalityPath[i,j] = @addConstraint(m, sum{t[totalPaths[i,j][a],totalPaths[i,j][a+1]], a=1:(length(totalPaths[i,j])-1)} - T[srcs[i],dsts[i]] >= - turnCost * totalNumExpensiveTurns[i,j] - delta * travelTimes[srcs[i],dsts[i]])
					end
				end
			end
		end
		# Equality constraints
		@addConstraint(m, equalityPath[i=1:numDataPoints], sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} - T[srcs[i],dsts[i]] == - turnCost * totalNumExpensiveTurns[i,1])
		toc()

		# Solve LP
		println("**** Solving LP ****")
		status = solve(m)
		if status == :Infeasible
			buildInternalModel(m)
			print_iis_gurobi(m)
			break
		end
		objective = getObjectiveValue(m)

		println("**** Setting up second LP ****")
		@addConstraint(m, fixObjective, sum{sqrt(numRides[i,j]/travelTimes[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]} <= 1.001*objective)
		@setObjective(m, Min, sum{delta2[i,j], i=nodes, j=out[i]})

		# Solve second LP
		println("**** Solving second LP ****")
		status = solve(m)

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
		if l <= max_rounds
			println("**** Computing shortest paths ****")
			@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
			@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
			compute_error(testingData, new_sp.traveltime)
		end
		l += 1
	end
	println("-- Saved outputs to directory Outputs/$(TESTDIR)/")
	return status, newTimes
end

manhattan = loadCityGraph()
if PREPROCESS
	@time outputPreprocessedConstraints(manhattan, "training", radius=RADIUS, numClusters=NUM_CLUSTERS, minRides=MIN_RIDES, sampleSize=SAMPLE_SIZE, overwrite = false, times=TIMES)
end
travel_times, num_rides = loadNewTravelTimeData(trainOrTest = "training", radius=RADIUS, times=TIMES, min_rides=MIN_RIDES, preprocess=PREPROCESS, num_clusters=NUM_CLUSTERS, sampleSize = SAMPLE_SIZE);
testing_data, numRides = loadNewTravelTimeData(trainOrTest="testing", radius = RADIUS, times = TIMES, preprocess = false, loadTestingMatrixDirectly = true);
@time status, new_times = fast_LP(manhattan, travel_times, num_rides, testing_data)
