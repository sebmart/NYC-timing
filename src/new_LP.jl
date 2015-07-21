# new_LP.jl
# New LP formulation for travel time estimation
# Authored by Arthur J Delarue on 7/2/15

MODEL = "strict"
MAX_ROUNDS = 10
MIN_RIDES = 1
RADIUS = 140
TIMES = "1214"

PREPROCESS = true
NUM_CLUSTERS = 50
SAMPLE_SIZE = 5000

TURN_COST = 10.
LAMBDA = 1e5
EPSILON_BOUND = [1/2000,1/5000,1/10000]

function new_LP(
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
	lambda::Float64=LAMBDA)

	graph = manhattan.network
	roadTimes = deepcopy(manhattan.roadTime)
	distances = manhattan.distances
	positions = manhattan.positions
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
	m = Model(solver=GurobiSolver(TimeLimit=10000, Method=2, Crossover=0, OutputFlag=1))
	m2 = Model(solver=GurobiSolver(TimeLimit=10000, Method=2, Crossover=0, BarHomogeneous=1,OutputFlag=1,InfUnbdInfo=1))

	# Update lower bounds by making very short roads instantaneous
	for i = nodes, j = out[i]
		if roadTimes[i,j] < 1
			roadTimes[i,j] = 0
		end
	end

	# Add one variable for each road, for each model
	@defVar(m, t[i=nodes,j=out[i]] >= roadTimes[i,j])
	@defVar(m2, t2[i=nodes,j=out[i]] >= roadTimes[i,j])

	# Find nonzero travel times
	pairs = [find(travel_times[i,:]) for i=nodes]
	# Epsilon is the variable used for the absolute value difference between T and \hat{T}
	@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)
	# T is \hat{T} in the notes
	@defVar(m, T[i=nodes, j=pairs[i]] >= 0)
	@addConstraint(m, TLowerBound[i=nodes,j=pairs[i]], T[i,j] - travelTimes[i,j] >= - epsilon[i,j])
	@addConstraint(m, TUpperBound[i=nodes,j=pairs[i]], T[i,j] - travelTimes[i,j] <= epsilon[i,j])
	# Add the relaxation variables for the relaxed model
	if model_type == "relaxed"
		@defVar(m, delta[i=nodes, j=pairs[i]] >= 0)
		@addConstraint(m, convergenceBound[i=nodes, j=pairs[i]], delta[i,j] <= 0)
	elseif model_type == "strict"
		# do nothing
	else
		# Give the user a warning and some time to rectify his mistake
		println("!!!! Invalid model type !!!!")
		pause(5)
	end

	# Define objective variables and constraints for m2
	@defVar(m2, delta2[i=nodes,j=out[i]] >= 0)
	@addConstraint(m2, objConstrLower[i=nodes,j=out[i]], -1 * t2[i,j]/distances[i,j] + 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t2[j,k], k = out[j]} + sum{1/distances[h,i] * t2[h,i], h=inn[i]}) <= delta2[i,j])
	@addConstraint(m2, objConstrUpper[i=nodes,j=out[i]], t2[i,j]/distances[i,j] - 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t2[j,k], k = out[j]} + sum{1/distances[h,i] * t2[h,i], h=inn[i]}) <= delta2[i,j])
	
	# Define objective function for programs
	if model_type == "strict"
		@setObjective(m, Min, sum{ sqrt(numRides[i,j]/travel_times[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]})
	else
		@setObjective(m, Min, sum{ sqrt(numRides[i,j]/travelTimes[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]} + lambda * sum{delta[i,j]/travelTimes[i,j], i=nodes, j=pairs[i]})
	end
	@setObjective(m2, Min, sum{delta2[i,j], i=nodes, j=out[i]})

	# Create handles for all constraints
	# This is necessary because of the complex way in which we generate constraints
	numConstraints = sum([length(pairs[i]) for i=nodes])
	if model_type == "relaxed"
		@defConstrRef higherRelaxed[1:numConstraints, 1:max_rounds]
	end
	@defConstrRef lower[1:numConstraints, 1:max_rounds]
	@defConstrRef higher[1:numConstraints, 1:max_rounds]
	@defConstrRef lower2[1:numConstraints, 1:max_rounds]
	@defConstrRef higher2[1:numConstraints, 1:max_rounds]
	# Create objects to keep track of constraints that are already in the set of constraints
	# Map from hashed paths to constraint indices
	hashedPathIndices = fill((Uint64 => Int64)[], numConstraints)
	# Set of hashed paths for quick lookup
	hashedPaths = fill(Set(Uint64[]), numConstraints)
	# Keep track of last path
	previousPaths = [0 for i = 1:numConstraints]
	# Same for second LP (only need the set for lookup/insertion, the map does not change)
	hashedPaths2 = fill(Set(Uint64[]), numConstraints)

	status = 0
	newTimes = roadTimes
	objective = 0

	# Compute shortest paths (with turn cost)
	println("**** Computing shortest paths ****")
	@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, positions, turn_cost=turnCost)
	old_nodes = getInverseMapping(new_nodes, nv(new_graph))
	@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)

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
	# Create array in which we will store all the expensive turns
	# Beware, dealing with this array has been the #1 source of bugs recently
	totalNumExpensiveTurns = Array(Int, (numConstraints, max_rounds))

	l = 1
	while l <= max_rounds
		actual_l = l
		# Update bound on relaxation variables
		if model_type == "relaxed"
			for i = nodes, j=pairs[i]
				chgConstrRHS(convergenceBound[i,j], (1/10)^(div(actual_l-1,length(EPSILON_BOUND)))*EPSILON_BOUND[((actual_l-1)%length(EPSILON_BOUND)) + 1])
			end
		end
		# setSolver(m, GurobiSolver(TimeLimit=10000, Method=2))
		println("###### ROUND $l ######")
		# Add path constraints
		println("**** Adding constraints ****")
		tic()
		paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, newTimes, new_edge_dists)
		for i=1:numConstraints
			if previousPaths[i] != 0
				chgConstrRHS(lower[i,previousPaths[i]], TMAX)
				if model_type == "relaxed"
					chgConstrRHS(higher[i,previousPaths[i]], -TMAX)
				end
			end
			# If path is already in model
			if hash(paths[i]) in hashedPaths[i]
				# Choose this path to be the good one
				index = hashedPathIndices[i][hash(paths[i])]
				chgConstrRHS(lower[i,index], - turnCost * numExpensiveTurns[i])
				if model_type == "relaxed"
					chgConstrRHS(higher[i,index], - turnCost * numExpensiveTurns[i])
				end
			# If this is the first path for this pair of nodes
			elseif length(hashedPaths[i]) == 0
				hashedPaths[i] = Set([hash(paths[i])])
				hashedPathIndices[i] = [hash(paths[i]) => 1]
				higher[i,1] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] >= - turnCost * numExpensiveTurns[i])
				lower[i,1] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] <= - turnCost * numExpensiveTurns[i])
				if model_type == "relaxed"
					higherRelaxed[i,1] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] + delta[srcs[i],dsts[i]] >= - turnCost * numExpensiveTurns[i])
				end
				totalNumExpensiveTurns[i,1] = numExpensiveTurns[i]
			# If this is a new path but not the first one
			else
				push!(hashedPaths[i], hash(paths[i]))
				len = length(hashedPaths[i])
				hashedPathIndices[i][hash(paths[i])] = len
				higher[i,len] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] >= - turnCost * numExpensiveTurns[i])
				lower[i,len] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] <= - turnCost * numExpensiveTurns[i])
				if model_type == "relaxed"
					higherRelaxed[i,len] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] + delta[srcs[i],dsts[i]] >= - turnCost * numExpensiveTurns[i])
				end
				totalNumExpensiveTurns[i,len] = numExpensiveTurns[i]
			end
		end
		toc()

		# Solve LP
		println("**** Solving LP ****")
		status = solve(m)

		println("**** Setting up second LP ****")
		# Get delta values
		result = getValue(T)
		TValues = zeros(length(nodes), length(nodes))
		for element in result
			TValues[element[1], element[2]] = element[3]
		end
		if model_type == "relaxed"
			deltaResult = getValue(delta)
			deltaValues = zeros(length(nodes), length(nodes))
			d = 0
			n = 0
			for element in deltaResult
				d += element[3]
				n += 1
				deltaValues[element[1], element[2]] = element[3]
			end
			# println("......... ", d/n, " ", n, " ", maximum(deltaValues))
		end

		println("**** Adding constraints ****")
		# Set up second LP, add constraints
		@time for i=1:numConstraints
			if previousPaths[i] != 0
				chgConstrRHS(lower2[i,previousPaths[i]], TMAX)
				for hashedPath in hashedPaths2[i]
					index = hashedPathIndices[i][hashedPath]
					if model_type == "relaxed"
						value = 0.99 * TValues[srcs[i],dsts[i]] - turnCost * totalNumExpensiveTurns[i,index] - deltaValues[srcs[i],dsts[i]]
					else #if split(model_type, "_")[end] == "strict"
						value = 0.99 * TValues[srcs[i],dsts[i]] - turnCost * totalNumExpensiveTurns[i,index]
					end
					chgConstrRHS(higher2[i,index], value)
				end
			end
			if hash(paths[i]) in hashedPaths2[i]
				index = hashedPathIndices[i][hash(paths[i])]
				chgConstrRHS(lower2[i,index], 1.01 * TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				chgConstrRHS(higher2[i,index], 0.99 * TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				previousPaths[i] = index
			elseif length(hashedPaths2[i]) == 0
				hashedPaths2[i] = Set([hash(paths[i])])
				higher2[i,1] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} >= 0.99 * TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				lower2[i,1] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} <= 1.01 * TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				previousPaths[i] = 1
			else
				push!(hashedPaths2[i], hash(paths[i]))
				len = length(hashedPaths2[i])
				higher2[i,len] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} >= 0.99 * TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				lower2[i,len] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} <= 1.01 * TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				previousPaths[i] = len
			end
		end

		# Solve second LP
		println("**** Solving second LP ****")
		status = solve(m2)

		# Debug if infeasible
		if status == :Infeasible || status == :InfeasibleOrUnbounded
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
			if abs(getObjectiveValue(m) - objective)/objective < 1e-10
				save("Outputs/$TESTDIR/end.jld", "num_iter", l)
				l = max_rounds
			else
				objective = getObjectiveValue(m)
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
	close(outputFile)
	println("-- Saved outputs to directory Outputs/$(TESTDIR)/")
	return status, newTimes
end

manhattan = loadCityGraph()
if PREPROCESS
	@time outputPreprocessedConstraints(manhattan, "training", radius=RADIUS, numClusters=NUM_CLUSTERS, minRides=MIN_RIDES, sampleSize=SAMPLE_SIZE, overwrite = false, times=TIMES)
end
travel_times, num_rides = loadNewTravelTimeData(trainOrTest = "training", radius=RADIUS, times=TIMES, min_rides=MIN_RIDES, preprocess=PREPROCESS, num_clusters=NUM_CLUSTERS, sampleSize = SAMPLE_SIZE);
testing_data, numRides = loadNewTravelTimeData(trainOrTest="testing", radius = RADIUS, times = TIMES, preprocess = false);
@time status, new_times = new_LP(manhattan, travel_times, num_rides, testing_data)
