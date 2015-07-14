# new_LP_test.jl
# New LP formulation for travel time estimation - test for metropolis
# Authored by Arthur J Delarue on 7/9/15

PROB = 0.7
MODEL = "metropolis_$(PROB)"
MAX_ROUNDS = 50
MIN_RIDES = 1

TURN_COST = 2.0

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
	turnCost::Float64=TURN_COST)
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
	m = Model(solver=GurobiSolver(TimeLimit=10000, Method=2, Crossover=0))
	m2 = Model(solver=GurobiSolver(TimeLimit=10000, Method=2, Crossover=0, BarHomogeneous=1,OutputFlag=1))

	# Add one variable for each road
	@defVar(m, t[i=nodes,j=out[i]] >= roadTimes[i,j])
	@defVar(m2, t2[i=nodes,j=out[i]] >= roadTimes[i,j])

	# Add regularization variables
	pairs = [find(travel_times[i,:]) for i=nodes]
	@defVar(m, epsilon[i=nodes,j=pairs[i]] >= 0)
	if split(model_type, "_")[end] == "strict"
		@defVar(m, T[i=nodes, j=pairs[i]])
		@addConstraint(m, TLowerBound[i=nodes,j=pairs[i]], T[i,j] - travelTimes[i,j] >= - epsilon[i,j])
		@addConstraint(m, TUpperBound[i=nodes,j=pairs[i]], T[i,j] - travelTimes[i,j] <= epsilon[i,j])
	end

	# Define objective variables and constraints for m2
	@defVar(m2, delta2[i=nodes,j=out[i]] >= 0)
	@addConstraint(m2, objConstrLower[i=nodes,j=out[i]], -1 * t2[i,j]/distances[i,j] + 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t2[j,k], k = out[j]} + sum{1/distances[h,i] * t2[h,i], h=inn[i]}) <= delta2[i,j])
	@addConstraint(m2, objConstrUpper[i=nodes,j=out[i]], t2[i,j]/distances[i,j] - 1/(length(inn[i]) + length(out[j])) * (sum{1/distances[j,k] * t2[j,k], k = out[j]} + sum{1/distances[h,i] * t2[h,i], h=inn[i]}) <= delta2[i,j])
	
	# Define objective function for programs
	@setObjective(m, Min, sum{ sqrt(numRides[i,j]/travel_times[i,j]) * epsilon[i,j], i=nodes, j=pairs[i]})
	@setObjective(m2, Min, sum{delta2[i,j], i=nodes, j=out[i]})

	# Create handles for all constraints
	numConstraints = sum([length(pairs[i]) for i=nodes])
	@defConstrRef lower[1:numConstraints, 1:max_rounds]
	@defConstrRef higher[1:numConstraints, 1:max_rounds]
	@defConstrRef lower2[1:numConstraints, 1:max_rounds]
	@defConstrRef higher2[1:numConstraints, 1:max_rounds]
	# Create arrays to keep track of constraints that are already in the set of constraints
	hashedPathList = fill(copy(Uint64[]), numConstraints)
	pathListToIteration = fill(copy(Int[]), numConstraints, max_rounds))

	status = 0
	newTimes = roadTimes
	objective = 0

	# Compute shortest paths (with turn cost)
	println("**** Computing shortest paths ****")
	@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, positions, turn_cost=turnCost)
	old_nodes = getInverseMapping(new_nodes, nv(new_graph))
	@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)

	# Run over all pairs of nodes that have data
	lowerBounds = Float64[]	
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
	totalNumExpensiveTurns = Array(Int, (numConstraints, max_rounds))

	l = 1
	while l <= max_rounds
		actual_l = l
		# setSolver(m, GurobiSolver(TimeLimit=10000, Method=2))
		println("###### ROUND $l ######")
		# Add path constraints
		println("**** Adding constraints ****")
		tic()
		paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, newTimes, new_edge_dists)
		totalNumExpensiveTurns[:,l] = numExpensiveTurns
		for i=1:numConstraints
			if split(model_type, "_")[end] == "strict"
				higher[i,l] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] >= - turnCost * numExpensiveTurns[i])
				lower[i,l] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - T[srcs[i],dsts[i]] <= - turnCost * numExpensiveTurns[i])
			else
				higher[i,l] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} + epsilon[srcs[i],dsts[i]] >= travelTimes[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				lower[i,l] = @addConstraint(m, sum{t[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} - epsilon[srcs[i],dsts[i]] <= travelTimes[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
			end
			if l > 1
				chgConstrRHS(lower[i,l-1], TMAX)
			end
		end
		toc()

		# Solve LP
		println("**** Solving LP ****")
		status = solve(m)

		println("**** Setting up second LP ****")
		# Get epsilon/delta values
		if split(model_type, "_")[end] == "strict"
			result = getValue(T)
			TValues = zeros(length(nodes), length(nodes))
			for element in result
				TValues[element[1], element[2]] = element[3]
			end
		else
			epsilonResult = getValue(epsilon)
			epsilonValues = zeros(length(nodes), length(nodes))
			for element in epsilonResult
				epsilonValues[element[1], element[2]] = 1.001 * element[3]
			end
		end

		println("**** Adding constraints ****")
		# Set up second LP, add constraints
		@time for i=1:numConstraints
			if split(model_type, "_")[end] == "strict"
				higher2[i,l] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} >= 0.999*TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				lower2[i,l] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} <= 1.001*TValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				if l > 1
					chgConstrRHS(lower2[i,l-1], TMAX)
					for j = 1:(l-1)
						value = 0.999*TValues[srcs[i],dsts[i]] - turnCost * totalNumExpensiveTurns[i,j]
						chgConstrRHS(higher2[i,j], value)
					end
				end
			else 
				higher2[i,l] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} >= travelTimes[srcs[i],dsts[i]] - epsilonValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				lower2[i,l] = @addConstraint(m2, sum{t2[paths[i][a],paths[i][a+1]], a=1:(length(paths[i])-1)} <= travelTimes[srcs[i],dsts[i]] + epsilonValues[srcs[i],dsts[i]] - turnCost * numExpensiveTurns[i])
				if l > 1
					chgConstrRHS(lower2[i,l-1], TMAX)
					for j = 1:(l-1)
						value = travelTimes[srcs[i],dsts[i]] - epsilonValues[srcs[i],dsts[i]] - turnCost * totalNumExpensiveTurns[i,j]
						chgConstrRHS(higher2[i,j], value)
					end
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

@time status, new_times = new_LP(graph, travel_times, num_rides, minRoadTime, distances, positions)
@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, new_times, positions, turn_cost=TURN_COST)
@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
print("")
# compute_error(real_times, new_sp.traveltime; num_nodes = nv(graph))