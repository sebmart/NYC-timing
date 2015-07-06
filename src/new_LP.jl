# new_LP.jl
# New LP formulation for travel time estimation
# Authored by Arthur J Delarue on 7/2/15

MODEL = "lin"
MAX_ROUNDS = 1
MIN_RIDES = 3
RADIUS = 40
TIMES = "1214"

PREPROCESS = true
NUM_CLUSTERS = 50
SAMPLE_SIZE = 50000

TURN_COST = 10.

LAMBDA = 1e-5

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
	turnCost::Float64=TURN_COST,
	lambda::Float64=LAMBDA)

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
	writeDataToFile("Outputs/$TESTDIR", model_type, max_rounds, min_rides, radius, preprocess, num_clusters, sample_size, turnCost, lambda)

	println("-- Saving outputs to directory Outputs/$(TESTDIR)/")
	# Create output csv file to save algorithm data
	outputFile = open("Outputs/$TESTDIR/algorithm_output.csv", "w")

	# Create JuMP model
	println("**** Creating LP instance ****")
	m = Model(solver=GurobiSolver(TimeLimit=200))#, Method=1))

	# Add one variable for each road
	@defVar(m, t[i=nodes,j=out[i]])

	# Add regularization variables
	pairs = [find(travel_times[i,:]) for i=nodes]
	@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)

	# Define objective function
	if model_type == "lin"
		@setObjective(m, Min, lambda * sum{ t[i,j]/distances[i,j], i=nodes, j=out[i] } + sum{ epsilon[i,j]/sqrt(travel_times[i,j]), i=nodes, j=pairs[i]})
	elseif model_type == "quad"
		@setObjective(m, Min, lambda * sum{ t[i,j] * t[i,j]/(distances[i,j] * distances[i,j]), i=nodes, j=out[i] } + sum{ epsilon[i,j]/sqrt(travelTimes[i,j]), i=nodes, j=pairs[i]})
	else
		println("!!!! INVALID MODEL !!!!")
		return 0,0
	end

	# Add bounds on variables
	for i in nodes, j in out[i]
		@addConstraint(m, t[i,j] >= roadTimes[i,j])
	end

	# Keep track of number of path constraints added at each round
	numPathConstraints = zeros(max_rounds)

	status = 0
	newTimes = roadTimes

	# Compute shortest paths (with turn cost)
	println("**** Computing shortest paths ****")
	@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
	old_nodes = getInverseMapping(new_nodes, nv(new_graph))
	@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)

	l = 1
	while l <= max_rounds
		setSolver(m, GurobiSolver(TimeLimit=200))#, Method=1))
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
		numPathConstraints[l] = length(lowerBounds)
		paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, newTimes, new_edge_dists)
		@addConstraint(m, higher[i=1:length(lowerBounds),j=l], sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} + turnCost * numExpensiveTurns[i] - lowerBounds[i] >= - epsilon[paths[i][1],paths[i][end]])
		@addConstraint(m, lower[i=1:length(lowerBounds),j=l], sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} + turnCost * numExpensiveTurns[i] - lowerBounds[i] <= epsilon[paths[i][1],paths[i][end]])
		if l > 1
			chgConstraintRHS(lower[i=1:numPathConstraints[l-1], j=(l-1)], TMAX)
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
			if l == max_rounds
				counter = 0
				epsi = getValue(epsilon)
				for element in epsi
					if element[3] != 0
						counter += 1
					end
				end
				println("VIOLATED CONSTRAINTS: ", counter)
				appendDataToFile("Outputs/$TESTDIR", counter)
			end
			st = getValue(t)
			newTimes = spzeros(length(nodes), length(nodes))
			for element in st
				newTimes[element[1], element[2]] = element[3]
			end
			# Save updated Manhattan road times to file
			saveRoadTimes(newTimes, "$TESTDIR/manhattan-times-$l")
			# Check if convergence criterion is satisfied
			if violated_constraints == 0
				println("###### OPTIMUM REACHED ######")
				break
			end
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
	@time outputPreprocessedConstraints(manhattan, radius=RADIUS, numClusters=NUM_CLUSTERS, minRides=MIN_RIDES, sampleSize=SAMPLE_SIZE, overwrite = false)
end
travel_times = loadTravelTimeData(radius=RADIUS, times=TIMES, min_rides=MIN_RIDES, preprocess=PREPROCESS, num_clusters=NUM_CLUSTERS, sampleSize = SAMPLE_SIZE);
@time status, new_times = new_LP(manhattan, travel_times)