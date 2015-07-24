# new_LP_test.jl
# New LP formulation for travel time estimation - test for metropolis
# Authored by Arthur J Delarue on 7/9/15

VARMAP = {
	:Basic => 0,
	:Superbasic => -3,
	:NonbasicAtUpper => -2,
	:NonbasicAtLower => -1
}

PROB = 0.9
MODEL = "metropolis_$(PROB)_strict"
MAX_ROUNDS = 5   #DON'T CHOOSE 30
MIN_RIDES = 1

TURN_COST = 2.0

LAMBDA = [1e3]#[1e9,1e3,2e3,5e3,1e4,2e4,5e4,1e5,2e5,5e5,1e6,2e6,5e6,1e25]
DELTA_BOUND = [0.6,0.5,0.4]

function new_LP(
	graph::SimpleGraph,
	travelTimes::Array{Float64,2},
	numRides::Array{Int,2},
	roadTimes::AbstractArray{Float64, 2},
	distances::AbstractArray{Float64, 2},
	positions::Array{Coordinates};
	model_type::String=MODEL,
	max_rounds::Int=MAX_ROUNDS,
	min_rides::Int=MIN_RIDES,
	turnCost::Float64=TURN_COST,
	delta_bound::Array{Float64}=DELTA_BOUND)
	"""
	Runs our iterative LP formulation on a Metropolis object
	"""

	roads = edges(graph)
	nodes = vertices(graph)
	out = [copy(out_neighbors(graph,i)) for i in nodes]
	inn = [copy(in_neighbors(graph,i)) for i in nodes]

	TMAX = 1000
	TESTDIR = "minr$(min_rides)_i$(max_rounds)_$(model_type)"

	# Create directory if necessary:
	if !isdir("Outputs/$TESTDIR")
		mkdir("Outputs/$TESTDIR")
	end

	# Create output csv file to save algorithm data
	println("-- Saving outputs to directory Outputs/$(TESTDIR)/")

	# Create JuMP model for LP and QP
	println("**** Creating LP instance ****")
	m = Model(solver=GurobiSolver(TimeLimit=10000, Method=2, Crossover=1, OutputFlag=1))
	m2 = Model(solver=GurobiSolver(TimeLimit=10000, Method=1, OutputFlag=1, InfUnbdInfo=1))

	# Add one variable for each road
	@defVar(m, t[i=nodes,j=out[i]] >= roadTimes[i,j])
	@defVar(m2, t2[i=nodes,j=out[i]] >= roadTimes[i,j])

	# Add decision variables for first LP
	pairs = [find(travel_times[i,:]) for i=nodes]
	@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)
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
	numDataPoints = sum([length(pairs[i]) for i=nodes])
	@defConstrRef path[1:numDataPoints, 1:max_rounds]
	@defConstrRef path2[1:numDataPoints, 1:max_rounds]

	# Create arrays to keep track of constraints that are already in the set of constraints
	# Map from hashed paths to constraint indices
	hashedPathIndices = fill((Uint64 => Int64)[], numDataPoints)
	# Set of hashed paths for quick lookup
	hashedPaths = fill(Set(Uint64[]), numDataPoints)
	# Keep track of last path
	previousPaths = [0 for i = 1:numDataPoints]
	# Same for second LP
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
	sizehint(srcs, 12000)
	sizehint(dsts, 12000)
	for i in nodes, j in nodes
		if travelTimes[i,j] > 0
			# Load info to reconstruct shortest paths in parallel
			push!(srcs, i)
			push!(dsts, j)
		end
	end
	totalNumExpensiveTurns = Array(Int, (numDataPoints, max_rounds))

	# Create warmstart array for LP2
	startingValues = zeros(2 * length(roads))
	startingBasis = zeros(Int, 2 * length(roads))

	l = 1
	while l <= max_rounds
		actual_l = l
		if split(model_type, "_")[end] == "relaxed"
			delta = (1/2)^(div(actual_l-1,length(delta_bound)))*delta_bound[((actual_l-1)%length(delta_bound)) + 1]
			# Avoid numerical issues
			if delta < 1e-4
				delta = 0
			end
		end
		# setSolver(m, GurobiSolver(TimeLimit=10000, Method=2))
		println("###### ROUND $l ######")

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
					if split(model_type, "_")[end] == "relaxed"
						pathLowerBounds[path[i,index].idx] = - turnCost * totalNumExpensiveTurns[i,index] - delta * travelTimes[srcs[i],dsts[i]]
					elseif split(model_type, "_")[end] == "strict"
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
		(cbasis, rbasis) = MathProgBase.getbasis(im)

		println("**** Setting up second LP ****")
		# Get delta values
		result = MathProgBase.getsolution(im)
		for i = nodes, j = out[i]
			startingValues[t2[i,j].col] = result[t[i,j].col]
			startingValues[delta2[i,j].col] = NaN
			startingBasis[t2[i,j].col] = VARMAP[cbasis[t[i,j].col]]
		end
		TValues = zeros(length(nodes), length(nodes))
		for i = nodes, j = pairs[i]
			TValues[i, j] = result[T[i,j].col]
		end
		if split(model_type, "_")[end] == "relaxed"	
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
		for i = 641:1280
			println(MathProgBase.getbasis(im2)[1][i], MathProgBase.getsolution(im2)[i])
		end
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
			saveRoadTimes(newTimes, "$TESTDIR/metropolis-times-$l")
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
			@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, positions, turn_cost=turnCost)
			@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
			println("**** Error calculation ****")
			printErrorStatsToFile("Outputs/$TESTDIR/errorstats-$(actual_l).txt", real_times, travel_times, new_sp.traveltime, numRides, meanRoadTime, newTimes, num_nodes = nv(graph))
		end
		l += 1
	end
	return status, newTimes
end

function compute_error(realData::AbstractArray{Float64, 2}, inputData::AbstractArray{Float64, 2}, algorithmData::AbstractArray{Float64}, numRides::AbstractArray{Int, 2}, realLinkTimes::AbstractArray{Float64,2}, newLinkTimes::AbstractArray{Float64,2}; num_nodes::Int = 192, out=Array{Int}[])
	"""
	Given a matrix of TABs from real Data, input Data and our algorithm, returns some statistical measures of error
	"""
	num_links = 0
	num_results = 0
	num_nonzero_results = 0
	average_squared_error = 0
	average_relative_error = 0
	average_bias = 0
	average_nonzero_input_error = 0
	average_nonzero_result_error = 0
	average_link_relative_error = 0
	average_link_bias = 0
	for i = 1:num_nodes, j=1:num_nodes
		if realData[i,j] > 0 && (length(out) == 0 || j in out[i])
			average_squared_error = (num_results * average_squared_error + (realData[i,j] - algorithmData[i,j]) ^ 2) / (num_results + 1)
			average_relative_error = (num_results * average_relative_error + abs(realData[i,j] - algorithmData[i,j])/realData[i,j])/(num_results + 1)
			average_bias = (num_results * average_bias + (algorithmData[i,j] - realData[i,j]))/(num_results + 1)
			num_results += 1
			if numRides[i,j] > 0
				average_nonzero_input_error = (num_nonzero_results * average_nonzero_input_error + abs(realData[i,j] - inputData[i,j])/realData[i,j])/(num_nonzero_results + 1)
				average_nonzero_result_error = (num_nonzero_results * average_nonzero_result_error + abs(realData[i,j] - algorithmData[i,j])/realData[i,j])/(num_nonzero_results + 1)
				num_nonzero_results += 1
			end			
		end
		if realLinkTimes[i,j] > 0
			average_link_relative_error = (num_links * average_link_relative_error + abs(newLinkTimes[i,j]-realLinkTimes[i,j])/realLinkTimes[i,j])/(num_links + 1)
			average_link_bias = (num_links * average_link_bias + newLinkTimes[i,j] - realLinkTimes[i,j])/(num_links + 1)
			num_links += 1
		end
	end
	println("------------------------------------------------")
	println("Average squared error: \t\t\t\t\t", average_squared_error)
	println("Average relative error: \t\t\t\t", average_relative_error)
	println("Average bias: \t\t\t\t\t\t", average_bias)
	println("Average original relative error on nonzero rides: \t", average_nonzero_input_error)
	println("Average result relative error on nonzero rides: \t", average_nonzero_result_error)
	println("Average relative error on link times: \t\t\t", average_link_relative_error)
	println("Average bias on link times:\t\t\t\t", average_link_bias)
	println("------------------------------------------------")
	return average_squared_error, average_relative_error, average_bias, average_nonzero_input_error, average_nonzero_result_error, average_link_relative_error, average_link_bias
end

function printErrorStatsToFile(fileName::String, realData::AbstractArray{Float64, 2}, inputData::AbstractArray{Float64, 2}, algorithmData::AbstractArray{Float64}, numRides::AbstractArray{Int, 2}, realLinkTimes::AbstractArray{Float64,2}, newLinkTimes::AbstractArray{Float64,2}; num_nodes::Int = 192, out=Array{Int}[])
	average_squared_error, average_relative_error, average_bias, average_nonzero_input_error, average_nonzero_result_error, average_link_relative_error, average_link_bias = compute_error(realData, inputData, algorithmData, numRides, realLinkTimes, newLinkTimes, num_nodes=num_nodes, out=out)
	"""
	Same as compute_error, except prints to file as well
	"""
	file = open(fileName, "w")
	write(file, string("Average squared error: \t\t\t\t\t", average_squared_error, "\n"))
	write(file, string("Average relative error: \t\t\t\t", average_relative_error, "\n"))
	write(file, string("Average bias: \t\t\t\t\t\t", average_bias, "\n"))
	write(file, string("Average original relative error on nonzero rides: \t", average_nonzero_input_error, "\n"))
	write(file, string("Average result relative error on nonzero rides: \t", average_nonzero_result_error, "\n"))
	write(file, string("Average relative error on link times: \t\t\t", average_link_relative_error, "\n"))
	write(file, string("Average bias on link times:\t\t\t\t", average_link_bias, "\n"))
	close(file)
end

# This code loads the data and calls the LP
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

travel_times = load("Inputs/input-travelTimes-$(PROB).jld", "travelTimes")
num_rides = load("Inputs/input-numRides-$(PROB).jld", "numRides")

function compute_cost_function(realData::AbstractArray{Float64,2}, inputData::AbstractArray{Float64,2}, numRides::AbstractArray{Int,2}; num_nodes = 192)
	cost = 0
	for i = 1:num_nodes, j=1:num_nodes
		if numRides[i,j] > 0
			cost += abs(realData[i,j] - inputData[i,j]) * sqrt(numRides[i,j]/inputData[i,j])
		end
	end
	return cost
end

@time status, new_times = new_LP(graph, travel_times, num_rides, minRoadTime, distances, positions)
@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, new_times, positions, turn_cost=TURN_COST)
@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
print("")
# compute_cost_function(real_times, travel_times, num_rides, num_nodes = nv(graph))