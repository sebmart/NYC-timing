# new_LP.jl
# New LP formulation for travel time estimation
# Authored by Arthur J Delarue on 7/2/15

MODEL = "quad"
MAX_ROUNDS = 1
MIN_RIDES = 3
RADIUS = 140
TIMES = "1214"

PREPROCESS = true
NUM_CLUSTERS = 50
SAMPLE_SIZE = 50000

TURN_COST = 10.

function new_LP(
	manhattan::Manhattan,
	travelTimes::Array{Float64,2};
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

	# Create JuMP model
	println("**** Creating LP instance ****")
	m = Model(solver=GurobiSolver(TimeLimit=10000, Method=2))

	# Add one variable for each road
	@defVar(m, t[i=nodes,j=out[i]] >= roadTimes[i,j])

	# Add regularization variables
	pairs = [find(travel_times[i,:]) for i=nodes]
	@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)

	# Define objective function
	@setObjective(m, Min, sum{ epsilon[i,j]/sqrt(travel_times[i,j]), i=nodes, j=pairs[i]})

	# Create handles for all constraints
	numConstraints = sum([length(pairs[i]) for i=nodes])
	@defConstrRef lower[1:numConstraints,1:max_rounds]
	@defConstrRef higher[1:numConstraints,1:max_rounds]

	status = 0
	newTimes = roadTimes

	# Compute shortest paths (with turn cost)
	println("**** Computing shortest paths ****")
	@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
	old_nodes = getInverseMapping(new_nodes, nv(new_graph))
	@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)

	l = 1
	while l <= max_rounds
		# setSolver(m, GurobiSolver(TimeLimit=10000, Method=2))
		println("###### ROUND $l ######")
		# Add path constraints
		println("**** Adding constraints ****")
		tic()
		# Run over all pairs of nodes that have data
		paths = Array{Float64}[]
		lowerBounds = Float64[]
		numExpensiveTurns = Int[]
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
		paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, newTimes, new_edge_dists)
		for i=1:numConstraints
			higher[i,l] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} + turnCost * numExpensiveTurns[i] + epsilon[paths[i][1],paths[i][end]] >= lowerBounds[i])
			lower[i,l] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} + turnCost * numExpensiveTurns[i] - epsilon[paths[i][1],paths[i][end]] <= lowerBounds[i])
		end
		if l > 1
			for i=1:numConstraints
				chgConstrRHS(lower[i, l-1], TMAX)
			end
		end
		toc()

		# Solve LP
		println("**** Solving LP ****")
		status = solve(m)

		# Debug if infeasible
		if status == :Infeasible
			buildInternalModel(m)
			print_iis_gurobi(m)
			break
		# Prepare output
		elseif status == :Optimal
			st = getValue(t)
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
		else
			epsilonValues = zeros(length(nodes), length(nodes))
			epsilonResult = getValue(epsilon)
			for element in epsilonResult
				epsilonValues[element[1], element[2]] = element[3]
			end
		end
		l += 1
	end

	# Enter second phase of LP, set new objective
	println("#### FINAL ROUND ####")
	if model_type == "lin"
		@setObjective(m, Min, sum{ t[i,j]/distances[i,j], i=nodes, j=out[i]})
	else
		@setObjective(m, Min, sum{ t[i,j] * t[i,j]/(distances[i,j] * distances[i,j]), i=nodes, j=out[i] })
	end

	# Reset solver
	setSolver(m, GurobiSolver(TimeLimit=1000))

	# Make sure epsilon are no longer decision variables
	@addConstraint(m, e[i=nodes, j=pairs[i]], epsilon[i,j] == epsilonValues[i,j])

	# Solve
	solve(m)

	# Get final edge times and save
	st = getValue(t)
	newTimes = spzeros(length(nodes), length(nodes))
	for element in st
		newTimes[element[1], element[2]] = element[3]
	end
	saveRoadTimes(newTimes, "$TESTDIR/manhattan-times-$(l+1)")

	close(outputFile)
	return status, newTimes
end

manhattan = loadCityGraph()
if PREPROCESS
	@time outputPreprocessedConstraints(manhattan, radius=RADIUS, numClusters=NUM_CLUSTERS, minRides=MIN_RIDES, sampleSize=SAMPLE_SIZE, overwrite = false)
end
travel_times = loadTravelTimeData(radius=RADIUS, times=TIMES, min_rides=MIN_RIDES, preprocess=PREPROCESS, num_clusters=NUM_CLUSTERS, sampleSize = SAMPLE_SIZE);
@time status, new_times = new_LP(manhattan, travel_times)