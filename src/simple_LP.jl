# simple_LP.jl
# Use our LP formulation with the constraint that all edges have the same velocity
# Authored by Arthur J Delarue on 7/31/15

function simple_LP(
	manhattan::Manhattan,
	travelTimes::Array{Float64,2},
	numRides::Array{Int,2},
	testingData::Array{Float64,2},
	directory::String;
	model_type::String="relaxed",
	max_rounds::Int=5,
	turnCost::Float64=10.0,
	delta_bound::Array{Float64}=[0.06,0.05,0.04],
	maxNumPathsPerOD::Int=3,
	computeFinalSP::Bool=false)

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

	TESTDIR = directory

	# Create directory if necessary:
	if !isdir("Outputs/$TESTDIR")
		mkdir("Outputs/$TESTDIR")
	end

	# Tell us where output is going
	println("-- Saving outputs to directory Outputs/$(TESTDIR)/")
	turnCostFile = open("Outputs/$(TESTDIR)/tc-simple.csv", "w")

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
	@time new_graph, new_edge_dists, new_nodes, new_edge_isExpensive = modifyGraphForDijkstra(graph, newTimes, positions, turn_cost=turnCost)
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
		@defVar(m, t[i=nodes,j=out[i]] >= 0.038 * manhattan.distances[i,j])

		# Add decision variable for turn cost and inverse velocity
		@defVar(m, tc >= 0)
		@defVar(m, w >= 0)
		@addConstraint(m, oneVelocity[i=nodes, j=out[i]], t[i,j] == w * distances[i,j])

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
		tic()
		paths, numExpensiveTurns = reconstructMultiplePathsWithExpensiveTurnsParallel(new_sp.previous, srcs, dsts, old_nodes, new_sp.real_destinations, new_edge_isExpensive)
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
						inequalityPath[i,j] = @addConstraint(m, sum{t[totalPaths[i,j+1][a],totalPaths[i,j+1][a+1]], a=1:(length(totalPaths[i,j+1])-1)} - sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} >= - tc * totalNumExpensiveTurns[i,j+1] + tc * totalNumExpensiveTurns[i,1])
					elseif model_type == "relaxed"
						inequalityPath[i,j] = @addConstraint(m, sum{t[totalPaths[i,j+1][a],totalPaths[i,j+1][a+1]], a=1:(length(totalPaths[i,j+1])-1)} - sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} >= - tc * totalNumExpensiveTurns[i,j+1] + tc * totalNumExpensiveTurns[i,1] - delta * travelTimes[srcs[i],dsts[i]])
					end
				end
			end
		end
		# Equality constraints (shortest path close to travel time data)
		@addConstraint(m, TLowerBound[i=1:numDataPoints], sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} + tc * totalNumExpensiveTurns[i,1] - travelTimes[srcs[i],dsts[i]] >= - epsilon[srcs[i],dsts[i]])
		@addConstraint(m, TUpperBound[i=1:numDataPoints], sum{t[totalPaths[i,1][a],totalPaths[i,1][a+1]], a=1:(length(totalPaths[i,1])-1)} + tc * totalNumExpensiveTurns[i,1] - travelTimes[srcs[i],dsts[i]] <= epsilon[srcs[i],dsts[i]])
		toc()

		# Solve LP
		println("**** Solving LP ****")
		status = solve(m)

		# Debug if infeasible
		if status == :Infeasible
			println("!!!! Computing IIS !!!!")
			print_iis_gurobi(m)
			break
		# Prepare output
		elseif status == :Optimal || status == :Suboptimal
			turnCost = getValue(tc)
			println(string("Left turn cost (s): ", turnCost))
			println(string("Average speed (kph): ", 3.6/getValue(w)))
			write(turnCostFile, string(l,",",turnCost, "\n"))
			result = getValue(t)
			newTimes = spzeros(length(nodes), length(nodes))
			for element in result
				newTimes[element[1], element[2]] = element[3]
			end
			# Save updated Manhattan road times to file
			saveRoadTimes(newTimes, "$TESTDIR/manhattan-simple-times-$l")
			objective = getObjectiveValue(m)
			if abs(objective - old_objective)/old_objective < 1e-3
				save("Outputs/$TESTDIR/simple-end.jld", "num_iter", l)
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
			@time new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=turnCost)
			@time new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
			compute_error(testingData, new_sp.traveltime)
		end
		l += 1
	end
	println("-- Saved outputs to directory Outputs/$(TESTDIR)/")
	close(turnCostFile)
	return status, newTimes, turnCost, totalPaths, totalNumExpensiveTurns
end