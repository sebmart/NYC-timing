# build_LP.jl
# Once Manhattan problem is loaded, builds instance of LP necessary
# Authored by Arthur J Delarue on 6/9/15

MODEL = "quad"
REGULARIZATION = true
MAX_ROUNDS = 1
MIN_RIDES = 3
RADIUS = 140
TIMES = "1214"

PREPROCESS = true
NUM_CLUSTERS = 50
SAMPLE_SIZE = 100000

TURN_COST = 10.

LAMBDA = 1e2
LAMBDA_VALUES = [0, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0]

function build_LP(
	manhattan::Manhattan,
	travelTimes::Array{Float64,2};
	model_type::String=MODEL,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	radius::Int=RADIUS,
	times::String=TIMES,
	regularization::Bool=REGULARIZATION,
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

	TESTDIR = "r$(radius)_minr$(min_rides)_i$(max_rounds)_wd_$(times)_$(model_type)"
	if regularization
		TESTDIR=string(TESTDIR,"_reg")
	end
	if preprocess
		TESTDIR=string(TESTDIR, "_clust$(num_clusters)_rides$(sample_size)")
	end

	# Create directory if necessary:
	if !isdir("Outputs/$TESTDIR")
		mkdir("Outputs/$TESTDIR")
	end
	# Add vital information to file
	writeDataToFile("Outputs/$TESTDIR", model_type, max_rounds, min_rides, radius, regularization, preprocess, num_clusters, sample_size, turnCost, lambda)

	println("-- Saving outputs to directory Outputs/$(TESTDIR)/")
	# Create output csv file to save algorithm data
	outputFile = open("Outputs/$TESTDIR/algorithm_output.csv", "w")

	# Create JuMP model
	println("**** Creating LP instance ****")
	m = Model(solver=GurobiSolver(TimeLimit=200, Method=1))

	# Add one variable for each road
	@defVar(m, t[i=nodes,j=out[i]])

	# Add regularization variables if necessary
	if regularization
		pairs = [find(travel_times[i,:]) for i=nodes]
		@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)
	end

	# Define objective function
	if model_type == "quad"  && !regularization
		@setObjective(m, Min, sum{ t[i,j] * t[i,j]/(distances[i,j] * distances[i,j]), i=nodes, j=out[i] })
	elseif model_type == "lin"
		@setObjective(m, Min, sum{ t[i,j], i=nodes, j=out[i] })
	elseif model_type == "quad" && regularization
		@setObjective(m, Min, sum{ t[i,j] * t[i,j]/(distances[i,j] * distances[i,j]), i=nodes, j=out[i] } + lambda * sum{ epsilon[i,j]/sqrt(travelTimes[i,j]), i=nodes, j=pairs[i]})
	else
		println("!!!! INVALID MODEL !!!!")
		return 0,0
	end

	# Add bounds on variables
	for i in nodes, j in out[i]
		@addConstraint(m, t[i,j] >= roadTimes[i,j])
	end

	status = 0
	newTimes = roadTimes

	# Compute shortest paths (with turn cost)
	println("**** Computing shortest paths ****")
	@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
	old_nodes = getInverseMapping(new_nodes, nv(new_graph))
	@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)

	l = 1
	while l <= max_rounds
		#setSolver(m, GurobiSolver(TimeLimit=200, Method=1))
		println("###### ROUND $l ######")
		# Add path constraints
		println("**** Adding constraints ****")
		tic()
		# Set up stat counters
		total_constraints = 0
		violated_constraints = 0
		sum_violated_constraints = 0
		# Run over all pairs of nodes that have data
		paths = Array{Float64}[]
		lowerBounds = Float64[]
		numExpensiveTurns = Int[]
		if l == 1
			srcs = Int[]
			dsts = Int[]
			sizehint(srcs, sample_size)
			sizehint(dsts, sample_size)
			sizehint(lowerBounds, sample_size)
			for i in nodes, j in nodes
				if travelTimes[i,j] > 0
					total_constraints += 1
					# See if constraint is violated and by how much
					if new_sp.traveltime[i,j] < travelTimes[i,j]
						violated_constraints += 1
						sum_violated_constraints += (travelTimes[i,j] - new_sp.traveltime[i,j])
						# Load info to reconstruct shortest paths in parallel
						push!(srcs, i)
						push!(dsts, j)
						push!(lowerBounds, travelTimes[i,j])
					end
				end
			end
			paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, newTimes, new_edge_dists)
		else
			for i in nodes, j in nodes
				if travelTimes[i,j] > 0
					total_constraints += 1
					# See if constraint is violated and by how much
					if new_sp.traveltime[i,j] < travelTimes[i,j]
						violated_constraints += 1
						sum_violated_constraints += (travelTimes[i,j] - new_sp.traveltime[i,j])
						# Reconstruct shortest path
						pathNodes, expensiveTurns = reconstructPathWithExpensiveTurns(new_sp.previous[i,:], i, j, old_nodes, new_sp.real_destinations[i,:], newTimes, new_edge_dists)
						push!(paths, pathNodes)
						push!(lowerBounds, travelTimes[i,j])
						push!(numExpensiveTurns, expensiveTurns)
					end
				end
			end
		end
		if regularization
			@addConstraint(m, c[i=1:violated_constraints], sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} >= lowerBounds[i] - epsilon[paths[i][1],paths[i][end]] - turnCost * numExpensiveTurns[i])
		else
			@addConstraint(m, c[i=1:violated_constraints], sum{t[paths[i][a],paths[i][a+1]],a=1:(length(paths[i])-1)} >= lowerBounds[i] - turnCost * numExpensiveTurns[i])
		end
		# Write stats to file
		if regularization
			write(outputFile,string(lambda, ","))
		end
		write(outputFile,string(l,",",total_constraints,",",violated_constraints,",",sum_violated_constraints/violated_constraints,"\n"))
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
			if regularization && l == max_rounds
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

function test_lambda_values(
	manhattan::Manhattan,
	travelTimes::Array{Float64,2},
	lambda_values::Array{Float64};
	appendToFile::Bool = false,
	model_type::String=MODEL,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	radius::Int=RADIUS,
	times::String=TIMES,
	regularization::Bool=REGULARIZATION,
	preprocess::Bool=PREPROCESS,
	num_clusters::Int=NUM_CLUSTERS,
	sample_size::Int=SAMPLE_SIZE,
	turnCost::Float64=TURN_COST)
	
	# Reconstruct directory name
	TESTDIR = "r$(radius)_minr$(min_rides)_i$(max_rounds)_wd_$(TIMES)_$(model_type)"
	if regularization
		TESTDIR=string(TESTDIR,"_reg")
	end
	if preprocess
		TESTDIR=string(TESTDIR, "_clust$(num_clusters)_rides$(sample_size)")
	end
	# Make new file that will only store relevant information
	outputFileName = "Outputs/lambda_values_$TESTDIR.csv"
	if appendToFile
		f = open(outputFileName, "a")
	else
		f = open(outputFileName, "w")
	end
	# Solve LP for different lambdas, and save results
	for lambda_value in lambda_values
		println("########################################")
		println("# Testing lambda = $lambda_value")
		println("########################################")
		status, new_times = build_LP(manhattan, travelTimes, model_type=model_type, max_rounds = max_rounds, min_rides = min_rides, radius = radius, times = times, regularization = regularization, preprocess = preprocess, num_clusters = num_clusters, sample_size = sample_size, turnCost = turnCost, lambda = lambda_value)
		write(f, string(lambda_value, ","))
		write(f, string(sum(new_times), ","))
		infoFile = open("Outputs/$TESTDIR/info.txt")
		lines = readlines(infoFile)
		write(f, string(strip(split(lines[end], ":")[2]), "\n"))
		close(infoFile)
	end
	close(f)
	return 0
end

manhattan = loadCityGraph()
if PREPROCESS
	@time outputPreprocessedConstraints(manhattan, radius=RADIUS, numClusters=NUM_CLUSTERS, minRides=MIN_RIDES, sampleSize=SAMPLE_SIZE, overwrite = false)
end
travel_times = loadTravelTimeData(radius=RADIUS, times=TIMES, min_rides=MIN_RIDES, preprocess=PREPROCESS, num_clusters=NUM_CLUSTERS, sampleSize = SAMPLE_SIZE);
# @time status, new_times = build_LP(manhattan, travel_times)
test_lambda_values(manhattan, travel_times, LAMBDA_VALUES, appendToFile = false)