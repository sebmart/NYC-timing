# new_LP.jl
# New LP formulation for travel time estimation
# Authored by Arthur J Delarue on 7/2/15

MODEL = "avgtr"
MAX_ROUNDS = 5
MIN_RIDES = 1
RADIUS = 140
TIMES = "1214"

PREPROCESS = true
NUM_CLUSTERS = 50
SAMPLE_SIZE = 50000

TURN_COST = 10.

function new_LP(
	manhattan::Manhattan,
	travelTimes::Array{Float64,2},
	numRides::Array{Int,2};
	model_type::String=MODEL,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	radius::Int=RADIUS,
	times::String=TIMES,
	preprocess::Bool=PREPROCESS,
	num_clusters::Int=NUM_CLUSTERS,
	sample_size::Int=SAMPLE_SIZE,
	turnCost::Float64=TURN_COST)

	graph = manhattan.network
	roadTimes = manhattan.roadTime
	distances = manhattan.distances
	roads = edges(graph)
	nodes = vertices(graph)
	out = [copy(out_neighbors(graph,i)) for i in nodes]
	inn = [copy(in_neighbors(graph,i)) for i in nodes]

	TMAX = 3600
	TESTDIR = "new_r$(radius)_minr$(min_rides)_i$(max_rounds)_wd_$(times)_$(model_type)"
	if preprocess
		TESTDIR=string(TESTDIR, "_clust$(num_clusters)_rides$(sample_size)")
	end

	# Create directory if necessary:
	if !isdir("Outputs/$TESTDIR")
		mkdir("Outputs/$TESTDIR")
	end
	# Add vital information to file
	writeDataToFile("Outputs/$TESTDIR", model_type, max_rounds, min_rides, radius, preprocess, num_clusters, sample_size, turnCost)

	# Create output csv file to save algorithm data
	println("-- Saving outputs to directory Outputs/$(TESTDIR)/")
	outputFile = open("Outputs/$TESTDIR/algorithm_output.csv", "w")

	# Create JuMP model for LP and QP
	println("**** Creating LP instance ****")
	m = Model(solver=GurobiSolver(TimeLimit=10000, Method=2, Crossover=0))
	m2 = Model(solver=GurobiSolver(TimeLimit=10000, Method=2, Crossover=0, BarHomogeneous=1, InfUnbdInfo=1))

	# Add one variable for each road
	@defVar(m, t[i=nodes,j=out[i]] >= roadTimes[i,j])
	@defVar(m2, t2[i=nodes,j=out[i]] >= roadTimes[i,j])

	# Add regularization variables
	pairs = [find(travel_times[i,:]) for i=nodes]
	@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)

	# Define objective variables and constraints for m2
	@defVar(m2, delta[i=nodes,j=out[i]] >= 0)
	@addConstraint(m2, objConstrLower[i=nodes,j=out[i]], -1 * t2[i,j]/distances[i,j] + 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t2[j,k], k = out[j]} + sum{1/distances[h,i] * t2[h,i], h=inn[i]}) <= delta[i,j])
	@addConstraint(m2, objConstrUpper[i=nodes,j=out[i]], t2[i,j]/distances[i,j] - 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t2[j,k], k = out[j]} + sum{1/distances[h,i] * t2[h,i], h=inn[i]}) <= delta[i,j])
	
	# Define objective function for programs
	@setObjective(m, Min, sum{ sqrt(numRides[i,j]/travel_times[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]})
	@setObjective(m2, Min, sum{delta[i,j], i=nodes, j=out[i]})

	# Create handles for all constraints
	numConstraints = sum([length(pairs[i]) for i=nodes])
	@defConstrRef lower[1:numConstraints,1:max_rounds]
	@defConstrRef higher[1:numConstraints,1:max_rounds]
	@defConstrRef lower2[1:numConstraints,1:max_rounds]
	@defConstrRef higher2[1:numConstraints,1:max_rounds]

	status = 0
	newTimes = roadTimes

	# Compute shortest paths (with turn cost)
	println("**** Computing shortest paths ****")
	@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
	old_nodes = getInverseMapping(new_nodes, nv(new_graph))
	@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)

	# Run over all pairs of nodes that have data
	lowerBounds = Float64[]	
	srcs = Int[]
	dsts = Int[]
	sizehint(srcs, sample_size)
	sizehint(dsts, sample_size)
	sizehint(lowerBounds, sample_size)
	for i in nodes, j in nodes
		if travelTimes[i,j] > 0
			# Load info to reconstruct shortest paths in parallel
			push!(srcs, i)
			push!(dsts, j)
			push!(lowerBounds, travelTimes[i,j])
		end
	end
	totalNumExpensiveTurns = Array(Int, (numConstraints, max_rounds))

	l = 1
	while l <= max_rounds
		# setSolver(m, GurobiSolver(TimeLimit=10000, Method=2))
		println("###### ROUND $l ######")
		# Add path constraints
		println("**** Adding constraints ****")
		tic()
		paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, newTimes, new_edge_dists)
		totalNumExpensiveTurns[:,l] = numExpensiveTurns
		for i=1:numConstraints
			higher[i,l] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} + epsilon[paths[i][1],paths[i][end]] >= lowerBounds[i] - turnCost * numExpensiveTurns[i])
			lower[i,l] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - epsilon[paths[i][1],paths[i][end]] <= lowerBounds[i] - turnCost * numExpensiveTurns[i])
			if l > 1
				chgConstrRHS(lower[i,l-1], TMAX)
			end
		end
		toc()

		# Solve LP
		println("**** Solving LP ****")
		status = solve(m)

		println("**** Setting up second LP ****")
		# Get epsilon values
		epsilonResult = getValue(epsilon)
		epsilonValues = zeros(length(nodes), length(nodes))
		for element in epsilonResult
			epsilonValues[element[1], element[2]] = element[3]
		end

		println("**** Adding constraints ****")
		# Set up second LP, add constraints
		@time for i=1:numConstraints
			higher2[i,l] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} >= lowerBounds[i] - epsilonValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
			lower2[i,l] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} <= lowerBounds[i] + epsilonValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
			if l > 1
				chgConstrRHS(lower2[i,l-1], TMAX)
				for j = 1:(l-1)
					value = lowerBounds[i] - epsilonValues[srcs[i],dsts[i]] - turnCost * totalNumExpensiveTurns[i,j]
					chgConstrRHS(higher2[i,j], value)
				end
			end
		
		end

		println("**** Solving second LP ****")
		# Solve second LP
		status = solve(m2)

		# Debug if infeasible
		if status == :Infeasible
			println("!!!! Diagnosis pending !!!!")
			buildInternalModel(m2)
			print_iis_gurobi(m2)
			break
		# Prepare output
		elseif status == :Optimal || status == :Suboptimal
			st = getValue(t2)
			newTimes = spzeros(length(nodes), length(nodes))
			for element in st
				newTimes[element[1], element[2]] = element[3]
			end
			# Save updated Manhattan road times to file
			saveRoadTimes(newTimes, "$TESTDIR/manhattan-times-$l")
		elseif status == :UserLimit
			println("!!!! User time limit exceeded !!!!")
			break
		end
		if l < max_rounds
			println("**** Computing shortest paths ****")
			@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
			@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
		end
		l += 1
	end
	close(outputFile)
	return status, newTimes
end

manhattan = loadCityGraph()
if PREPROCESS
	@time outputPreprocessedConstraints(manhattan, "training", radius=RADIUS, numClusters=NUM_CLUSTERS, minRides=MIN_RIDES, sampleSize=SAMPLE_SIZE, overwrite = false, times="0709")
end
travel_times, num_rides = loadNewTravelTimeData(trainOrTest = "training", radius=RADIUS, times="0709", min_rides=MIN_RIDES, preprocess=PREPROCESS, num_clusters=NUM_CLUSTERS, sampleSize = SAMPLE_SIZE);
@time status, new_times = new_LP(manhattan, travel_times, num_rides)