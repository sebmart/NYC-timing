# new_LP.jl
# New LP formulation for travel time estimation
# Authored by Arthur J Delarue on 7/2/15

MODEL = "relaxed"
MAX_ROUNDS = 10
MIN_RIDES = 1
RADIUS = 140
TIMES = "1214"

PREPROCESS = true
NUM_CLUSTERS = 50
SAMPLE_SIZE = 5000

TURN_COST = 10.

DELTA_BOUND = [0.06,0.05,0.04]

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
	delta_bound::Array{Float64}=DELTA_BOUND)

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
	m = Model(solver=GurobiSolver(TimeLimit=10000, Method=2, Crossover=1, OutputFlag=1))
	m2 = Model(solver=GurobiSolver(TimeLimit=10000, Method=1, InfUnbdInfo=1))

	# Update lower bounds by making very short roads instantaneous
	# for i = nodes, j = out[i]
	# 	if roadTimes[i,j] < 1
	# 		roadTimes[i,j] = 0
	# 	end
	# end

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

	# Define objective variables and constraints for m2
	@defVar(m2, delta2[i=nodes,j=out[i]] >= 0)
	@addConstraint(m2, objConstrLower[i=nodes,j=out[i]], -1 * t2[i,j]/distances[i,j] + 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t2[j,k], k = out[j]} + sum{1/distances[h,i] * t2[h,i], h=inn[i]}) <= delta2[i,j])
	@addConstraint(m2, objConstrUpper[i=nodes,j=out[i]], t2[i,j]/distances[i,j] - 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t2[j,k], k = out[j]} + sum{1/distances[h,i] * t2[h,i], h=inn[i]}) <= delta2[i,j])
	
	# Define objective function for programs
	@setObjective(m, Min, sum{ sqrt(numRides[i,j]/travel_times[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]})
	@setObjective(m2, Min, sum{delta2[i,j], i=nodes, j=out[i]})

	# Create handles for all constraints
	# This is necessary because of the complex way in which we generate constraints
	numDataPoints = sum([length(pairs[i]) for i=nodes])
	@defConstrRef path[1:numDataPoints, 1:max_rounds]
	@defConstrRef path2[1:numDataPoints, 1:max_rounds]

	# Create objects to keep track of constraints that are already in the set of constraints
	# Map from hashed paths to constraint indices
	hashedPathIndices = fill((Uint64 => Int64)[], numDataPoints)
	# Set of hashed paths for quick lookup
	hashedPaths = fill(Set(Uint64[]), numDataPoints)
	# Keep track of last path
	previousPaths = [0 for i = 1:numDataPoints]
	# Same for second LP (only need the set for lookup/insertion, the map does not change)
	hashedPaths2 = fill(Set(Uint64[]), numDataPoints)

	status = 0
	newTimes = roadTimes
	objective = 0
	totalPathConstraints = 0

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
	totalNumExpensiveTurns = Array(Int, (numDataPoints, max_rounds))

	# Create warmstart array for LP2
	startingValues = zeros(2 * length(roads))

	l = 1
	while l <= max_rounds
		println("###### ROUND $l ######")
		actual_l = l
		# Update bound on relaxation variables
		if model_type == "relaxed"
			delta = (1/2)^(div(actual_l-1,length(delta_bound)))*delta_bound[((actual_l-1)%length(delta_bound)) + 1]
			# Avoid numerical issues
			if delta < 1e-4
				delta = 0
			end
		end
		setSolver(m2, GurobiSolver(TimeLimit=10000, Method=1))

		# Create upper and lower bound arrays for MathProgBase manual constraint management
		pathLowerBounds = zeros(totalPathConstraints + 2 * numDataPoints)
		pathUpperBounds = zeros(totalPathConstraints + 2 * numDataPoints)
		# Fill in upper and lower bounds for non-path constraints
		for i = 1:numDataPoints
			pathLowerBounds[i] = travelTimes[srcs[i],dsts[i]]
			pathLowerBounds[numDataPoints + i] = -Inf
			pathUpperBounds[i] = Inf
			pathUpperBounds[numDataPoints + i] = travelTimes[srcs[i],dsts[i]]
		end

		# Add path constraints
		println("**** Adding constraints ****")
		tic()
		paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, newTimes, new_edge_dists)
		for i=1:numDataPoints
			# Update old paths
			if previousPaths[i] != 0
				for hashedPath in hashedPaths[i]
					index = hashedPathIndices[i][hashedPath]
					if model_type == "relaxed"
						pathLowerBounds[path[i,index].idx] = - turnCost * totalNumExpensiveTurns[i,index] - delta * travelTimes[srcs[i],dsts[i]]
					elseif model_type == "strict"
						pathLowerBounds[path[i,index].idx] = - turnCost * totalNumExpensiveTurns[i,index]
					end
					pathUpperBounds[path[i,index].idx] = Inf
				end
			end
			# If path is already in model
			if hash(paths[i]) in hashedPaths[i]
				# Choose this path to be the good one
				index = hashedPathIndices[i][hash(paths[i])]
				pathLowerBounds[path[i,index].idx] = - turnCost * numExpensiveTurns[i]
				pathUpperBounds[path[i,index].idx] = - turnCost * numExpensiveTurns[i]
			# If this is the first path for this pair of nodes
			elseif length(hashedPaths[i]) == 0
				hashedPaths[i] = Set([hash(paths[i])])
				hashedPathIndices[i] = [hash(paths[i]) => 1]
				path[i,1] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] == - turnCost * numExpensiveTurns[i])
				push!(pathLowerBounds, - turnCost * numExpensiveTurns[i])
				push!(pathUpperBounds, - turnCost * numExpensiveTurns[i])
				totalNumExpensiveTurns[i,1] = numExpensiveTurns[i]
				totalPathConstraints += 1
			# If this is a new path but not the first one
			else
				push!(hashedPaths[i], hash(paths[i]))
				len = length(hashedPaths[i])
				hashedPathIndices[i][hash(paths[i])] = len
				path[i,len] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] == - turnCost * numExpensiveTurns[i])
				push!(pathLowerBounds, - turnCost * numExpensiveTurns[i])
				push!(pathUpperBounds, - turnCost * numExpensiveTurns[i])
				totalNumExpensiveTurns[i,len] = numExpensiveTurns[i]
				totalPathConstraints += 1
			end
		end
		toc()

		# Solve LP
		println("**** Solving LP ****")
		buildInternalModel(m)
		im = getInternalModel(m)
		MathProgBase.setconstrLB!(im, pathLowerBounds)
		MathProgBase.setconstrUB!(im, pathUpperBounds)
		MathProgBase.updatemodel!(im)
		MathProgBase.optimize!(im)
		status = MathProgBase.status(im)

		println("**** Setting up second LP ****")
		# Get delta and T values + edge weights
		result = MathProgBase.getsolution(im)
		for i = nodes, j = out[i]
			startingValues[t2[i,j].col] = result[t[i,j].col]
			startingValues[delta2[i,j].col] = NaN
		end
		TValues = zeros(length(nodes), length(nodes))
		for i = nodes, j = pairs[i]
			TValues[i, j] = result[T[i,j].col]
		end
		if model_type == "relaxed"	
			# Find delta values from how much constraints were actually violated
			deltaValues = zeros(length(nodes), length(nodes))
			deltaPercentage = zeros(length(nodes), length(nodes))
			# Get values of vector Ax
			constraintValues = MathProgBase.getconstrsolution(im)
			for i = 1:numDataPoints
				difference = 0.0
				# Look at all paths (except equality)
				for hashedPath in hashedPaths[i]
					if hashedPath != hash(paths[i])
						# Find value of violation. deltaValue is maximum violation
						index = hashedPathIndices[i][hashedPath]
						newDiff = - constraintValues[path[i,index].idx] - totalNumExpensiveTurns[i,index] * turnCost
						if newDiff > difference
							difference = newDiff
						end
					end
				end
				deltaValues[srcs[i], dsts[i]] = difference
				deltaPercentage[srcs[i], dsts[i]] = difference/travelTimes[srcs[i],dsts[i]]
			end
			println("Max delta value: ", maximum(deltaPercentage))
		end

		# Create upper and lower bound arrays for MathProgBase manual constraint management
		pathLowerBounds2 = zeros(totalPathConstraints + 2 * length(roads))
		pathUpperBounds2 = zeros(totalPathConstraints + 2 * length(roads))
		for i = 1:2*length(roads)
			pathLowerBounds2[i] = -Inf
			pathUpperBounds2[i] = 0.0
		end

		println("**** Adding constraints ****")
		# Set up second LP, add constraints
		@time for i=1:numDataPoints
			if previousPaths[i] != 0
				for hashedPath in hashedPaths2[i]
					index = hashedPathIndices[i][hashedPath]
					if split(model_type, "_")[end] == "relaxed"
						value = TValues[srcs[i],dsts[i]] - turnCost * totalNumExpensiveTurns[i,index] - deltaValues[srcs[i],dsts[i]]
					else #if split(model_type, "_")[end] == "strict"
						value = TValues[srcs[i],dsts[i]] - turnCost * totalNumExpensiveTurns[i,index]
					end
					pathLowerBounds2[path2[i,index].idx] = value
					pathUpperBounds2[path2[i,index].idx] = Inf
				end
			end
			if hash(paths[i]) in hashedPaths2[i]
				index = hashedPathIndices[i][hash(paths[i])]
				pathLowerBounds2[path2[i,index].idx] = TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i]
				pathUpperBounds2[path2[i,index].idx] = TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i]
				previousPaths[i] = index
			elseif length(hashedPaths2[i]) == 0
				hashedPaths2[i] = Set([hash(paths[i])])
				path2[i,1] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} == TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				pathLowerBounds2[path2[i,1].idx] = TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i]
				pathUpperBounds2[path2[i,1].idx] = TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i]
				previousPaths[i] = 1
			else
				push!(hashedPaths2[i], hash(paths[i]))
				len = length(hashedPaths2[i])
				path2[i,len] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} == TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				pathLowerBounds2[path2[i,len].idx] = TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i]
				pathUpperBounds2[path2[i,len].idx] = TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i]
				previousPaths[i] = len
			end
		end

		# Solve second LP
		println("**** Solving second LP ****")
		buildInternalModel(m2)
		im2 = getInternalModel(m2)
		MathProgBase.setwarmstart!(im2, startingValues)
		MathProgBase.setconstrLB!(im2, pathLowerBounds2)
		MathProgBase.setconstrUB!(im2, pathUpperBounds2)
		MathProgBase.updatemodel!(im2)
		MathProgBase.optimize!(im2)
		status = MathProgBase.status(im2)

		# Debug if infeasible
		if status == :Infeasible
			println("!!!!!!!!!!!!!!!!!!!!!!!")
			println("!!!! Computing IIS !!!!")
			println("!!!!!!!!!!!!!!!!!!!!!!!")
			print_iis_gurobi(m2, im2)
			break
		# Prepare output
		elseif status == :Optimal || status == :Suboptimal
			result2 = MathProgBase.getsolution(im2)
			newTimes = spzeros(length(nodes), length(nodes))
			for i = nodes, j=out[i]
				newTimes[i, j] = result2[t2[i,j].col]
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
testing_data, numRides = loadNewTravelTimeData(trainOrTest="testing", radius = RADIUS, times = TIMES, preprocess = false, loadTestingMatrixDirectly = true);
@time status, new_times = new_LP(manhattan, travel_times, num_rides, testing_data)
