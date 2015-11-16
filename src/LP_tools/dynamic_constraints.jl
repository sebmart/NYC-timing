# dynamic_constraints.jl
# Dynamic constraint generation helper functions
# Authored by Arthur J Delarue on 9/27/15

MANHATTAN_NODES = 6134

function chooseStartingConstraints(
	travelTimes::Array{Float64, 2},
	numRides::Array{Int, 2},
	pairs::Array{Array{Int,1},1};
	sample_size::Int = 500,
	num_nodes::Int = MANHATTAN_NODES)
	"""
	Chooses constraints that algorithm will start with. Takes in matrix array of travel times and number of rides (and starting sample size) and returns arrays of travel times/number of rides of same decisions but only relevant info, plus 1D arrays with sources and destinations.
	"""
	# Randomly pick starting constraints
	newTravelTimes, newNumRides = chooseConstraints(travelTimes, numRides, sample_size = sample_size, num_nodes=num_nodes)
	# Create vectors for sources and destinations
	srcs = Int[]
	dsts = Int[]
	sizehint!(srcs, sample_size)
	sizehint!(dsts, sample_size)
	for i in 1:size(travelTimes)[1], j in 1:size(travelTimes)[2]
		if newTravelTimes[i,j] > 0
			# Load sources and destinations: will be useful when adding constraints
			push!(srcs, i)
			push!(dsts, j)
			push!(pairs[i], j)
		end
	end
	return srcs, dsts, newTravelTimes, newNumRides, pairs
end

function updateConstraints(
	travelTimes::Array{Float64,2},
	numRides::Array{Int,2},
	calculatedTravelTimes::Array{Float64,2},
	totalPaths::Array{Any, 1},
	totalNumExpensiveTurns::Array{Any, 1},
	numPaths::Array{Int, 1},
	srcs::Array{Int, 1},
	dsts::Array{Int, 1},
	pairs::Array{Array{Int, 1}, 1};
	numNodePairsToAdd::Int = 0,
	numNodePairsToRemove::Int = 0,
	addOnlyIfAbovePercentile::Float64 = 0.5,
	removeOnlyIfBelowPercentile::Float64 = 0.5)

	"""
	Add and remove node pairs based on selection criteria.
	"""

	# Create new return structures
	newSrcs = deepcopy(srcs)
	newDsts = deepcopy(dsts)
	newPairs = deepcopy(pairs)
	newTotalPaths = deepcopy(totalPaths)
	newTotalNumExpensiveTurns = deepcopy(totalNumExpensiveTurns)
	newNumPaths = deepcopy(numPaths)

	# Create lookup set for constraints already under consideration
	nodePairsAlreadyIn = Set([(srcs[i], dsts[i]) for i = 1:length(srcs)])
	# Select new constraints - first, put everything in vector form to prepare for sorting
	indicesVector = Tuple{Int,Int}[]
	sizehint!(indicesVector, 12000000)
	for i = 1:size(travelTimes)[1], j = 1:size(travelTimes)[1]
		if travelTimes[i,j] > 0
			push!(indicesVector, (i,j))
		end
	end
	# then, compute distance between calculated times and data for all pairs
	newErrorVector = zeros(length(indicesVector))
	for i = 1:length(indicesVector)
		newErrorVector[i] = sqrt(numRides[indicesVector[i][1], indicesVector[i][2]]/travelTimes[indicesVector[i][1], indicesVector[i][2]]) * abs(travelTimes[indicesVector[i][1], indicesVector[i][2]] - calculatedTravelTimes[indicesVector[i][1], indicesVector[i][2]])
	end
	# compute error for nodepairs already in LP
	errorVector = zeros(length(srcs))
	for i = 1:length(srcs)
		errorVector[i] = sqrt(numRides[srcs[i], dsts[i]]/travelTimes[srcs[i], dsts[i]]) * abs(travelTimes[srcs[i], dsts[i]] - calculatedTravelTimes[srcs[i], dsts[i]])
	end
	# sort the ones already in the data by error
	p2 = sortperm(errorVector)
	# sort everything by error
	p = sortperm(newErrorVector)
	# Find ranges where it is acceptable to add/remove constraints
	minIndexToAdd = round(Int, addOnlyIfAbovePercentile * length(indicesVector))
	maxIndexToRemove = round(Int, removeOnlyIfBelowPercentile * length(newSrcs))
	# Now, remove best node pairs - create vector of indices to keep, then splice out the indices we don't want to keep
	indicesToKeep = collect(1:length(newSrcs))
	numNodePairsRemoved = 0
	while numNodePairsRemoved < min(numNodePairsToRemove, maxIndexToRemove)
		i = rand(1:maxIndexToRemove)
		if (newSrcs[p2[i]], newDsts[p2[i]]) in nodePairsAlreadyIn
			idx = findfirst(indicesToKeep, p2[i])
			splice!(indicesToKeep, idx)
			idx = findfirst(newPairs[newSrcs[p2[i]]], newDsts[p2[i]])
			splice!(newPairs[newSrcs[p2[i]]], idx)
			delete!(nodePairsAlreadyIn, (newSrcs[p2[i]], newDsts[p2[i]]))
			numNodePairsRemoved += 1
		end
	end
	# Since we spliced out the unwanted indices, this performs the necessary reduction
	newSrcs = newSrcs[indicesToKeep]
	newDsts = newDsts[indicesToKeep]
	newTotalNumExpensiveTurns = newTotalNumExpensiveTurns[indicesToKeep]
	newNumPaths = newNumPaths[indicesToKeep]
	newTotalPaths = newTotalPaths[indicesToKeep]
	# Add poor node pairs - this is easier, just push the new nodepairs to the end :)
	numNodePairsAdded = 0
	while numNodePairsAdded < numNodePairsToAdd && length(nodePairsAlreadyIn) < (length(indicesVector) - minIndexToAdd)
		i = rand(minIndexToAdd:length(indicesVector))
		if !(indicesVector[p[i]] in nodePairsAlreadyIn)
			push!(newSrcs, indicesVector[p[i]][1])
			push!(newDsts, indicesVector[p[i]][2])
			push!(newTotalPaths, Array{Int}[])
			push!(newTotalNumExpensiveTurns, Int[])
			push!(newNumPaths, 0)
			push!(newPairs[indicesVector[p[i]][1]], indicesVector[p[i]][2])
			sort!(newPairs[indicesVector[p[i]][1]])
			numNodePairsAdded += 1
			push!(nodePairsAlreadyIn, indicesVector[p[i]])
		end
	end
	return newSrcs, newDsts, newTotalPaths, newTotalNumExpensiveTurns, newNumPaths, newPairs
end

function updatePaths(
	paths::Union{Array{Any, 1}, Array{Array{Int,1},1}},
	numExpensiveTurns::Array{Int,1},
	totalPaths::Union{Array{Any, 1}, Array{Array{Int,1},1}},
	totalNumExpensiveTurns::Array{Any,1},
	numPaths::Array{Int, 1};
	maxNumPathsPerOD::Int=3,
	times::AbstractArray{Float64,2}=zeros(1,1),
	turnCost::Float64=0.0,
	travelTimes::AbstractArray{Float64,2}=zeros(1,1),
	numRides::AbstractArray{Int,2}=zeros(Int,(1,1)),
	srcs::Array{Int,1}=zeros(1),
	dsts::Array{Int,1}=zeros(1),
	dynamicConstraints::Bool=true,
	globalConstraintUpdate::Bool=true
	)
	"""
	Update path array so as not to exceed max number of paths.
	"""
	newTotalPaths = [totalPaths[i] for i = 1:length(totalPaths)]
	newTotalNumExpensiveTurns = [totalNumExpensiveTurns[i] for i = 1:length(totalNumExpensiveTurns)]
	totalPaths = newTotalPaths
	totalNumExpensiveTurns = newTotalNumExpensiveTurns
	if dynamicConstraints && globalConstraintUpdate
		# Deal with single path per O,D case separately because of the bugs it's been causing
		if maxNumPathsPerOD > 1
			for i=1:length(paths)
				index = findfirst(totalPaths[i], paths[i])
				# If path already in paths
				if index != 0
					totalPaths[i][1], totalPaths[i][index] = totalPaths[i][index], totalPaths[i][1]
					totalNumExpensiveTurns[i][1], totalNumExpensiveTurns[i][index] = totalNumExpensiveTurns[i][index], totalNumExpensiveTurns[i][1]			# New path is equality constraint
				# If we still have space to add the path
				else
					numPaths[i] += 1
					if numPaths[i] == 1
						push!(totalPaths[i], paths[i])
						push!(totalNumExpensiveTurns[i], numExpensiveTurns[i])
					else
						push!(totalPaths[i], paths[i])
						push!(totalNumExpensiveTurns[i], numExpensiveTurns[i])
						assert(numPaths[i] == length(totalPaths[i]))
						totalPaths[i][numPaths[i]], totalPaths[i][1] = totalPaths[i][1], totalPaths[i][numPaths[i]]
						totalNumExpensiveTurns[i][numPaths[i]], totalNumExpensiveTurns[i][1] = totalNumExpensiveTurns[i][1], totalNumExpensiveTurns[i][numPaths[i]]
					end
				end
			end
			# Check if paths need to be removed, and if so, remove them
			nPaths = sum(numPaths)
			# Check if we have too many paths
			if nPaths > maxNumPathsPerOD * length(totalPaths)
				# Put all the paths that could be removed (i.e. not equality paths)
				# into a big array (indicesVector)
				# so that we will be able to find them after we sort
				indicesVector = Tuple{Int,Int,Int,Int}[]
				sizehint!(indicesVector, nPaths)
				for i = 1:length(totalPaths), j=1:numPaths[i]
					if j > 1
						push!(indicesVector, (i,j,srcs[i],dsts[i]))
					end
				end
				errorVector = zeros(length(indicesVector))
				# Check if there are any non-equality paths
				if length(indicesVector) > 0
					# Compute error on each path (difference between it and the equality path,
					# 						weighted by path time)
					for i=1:length(indicesVector)
						errorVector[i] = abs(sum([times[totalPaths[indicesVector[i][1]][indicesVector[i][2]][a], totalPaths[indicesVector[i][1]][indicesVector[i][2]][a+1]] for a = 1:(length(totalPaths[indicesVector[i][1]][indicesVector[i][2]])-1)]) + totalNumExpensiveTurns[indicesVector[i][1]][indicesVector[i][2]] * turnCost - sum([times[totalPaths[indicesVector[i][1]][1][a], totalPaths[indicesVector[i][1]][1][a+1]] for a = 1:(length(totalPaths[indicesVector[i][1]][1])-1)]) - totalNumExpensiveTurns[indicesVector[i][1]][1] * turnCost)/(travelTimes[indicesVector[i][3], indicesVector[i][4]])
					end
					# Sort paths by "badness"
					p = sortperm(errorVector)
					# Figure out how many paths need to be removed
					numPathsToRemove = nPaths - maxNumPathsPerOD * length(totalPaths)
					# Initialize an array of indices of the paths that need to be kept (so far all paths are in there)
					pathsToKeep = [collect(1:numPaths[i]) for i=1:length(numPaths)]
					for i = 1:numPathsToRemove
						# Remove appropriate path from index list
						idx = findfirst(pathsToKeep[indicesVector[p[end+1-i]][1]], indicesVector[p[end+1-i]][2])
						splice!(pathsToKeep[indicesVector[p[end+1-i]][1]], idx)
					end
					# Rebuild path arrays using pathsToKeep to filter out removed paths
					totalPaths = [[totalPaths[i][j] for j in pathsToKeep[i]] for i = 1:length(totalPaths)]
					totalNumExpensiveTurns = [[totalNumExpensiveTurns[i][j] for j in pathsToKeep[i]] for i = 1:length(totalNumExpensiveTurns)]
					numPaths = [length(pathsToKeep[i]) for i = 1:length(numPaths)]
				end
			end
		else
			# Recreate path and turn arrays from scratch (since we will overwrite all their contents anyways)
			numDataPoints = length(totalPaths)
			assert(numDataPoints == length(paths))
			totalPaths = Array(Any, numDataPoints)
			totalNumExpensiveTurns = Array(Any, numDataPoints)
			for i = 1:numDataPoints
				totalPaths[i] = Array{Int}[]
				totalNumExpensiveTurns[i] = Int[]
			end
			# If one path per constraint, just update the path... easy-peasy
			for i=1:length(paths)
				push!(totalPaths[i], paths[i])
				push!(totalNumExpensiveTurns[i], numExpensiveTurns[i])
			end
		end
	# Non-dynamic constraint case
	else
		for i=1:length(paths)
			index = findfirst(totalPaths[i], paths[i])
			# If path already in paths
			if index != 0
				totalPaths[i][1], totalPaths[i][index] = totalPaths[i][index], totalPaths[i][1]
				totalNumExpensiveTurns[i][1], totalNumExpensiveTurns[i][index] = totalNumExpensiveTurns[i][index], totalNumExpensiveTurns[i][1]			# New path is equality constraint
			# If we still have space to add the path
			elseif numPaths[i] < maxNumPathsPerOD
				numPaths[i] += 1
				if numPaths[i] == 1
					push!(totalPaths[i], paths[i])
					push!(totalNumExpensiveTurns[i], numExpensiveTurns[i])
				else
					push!(totalPaths[i], paths[i])
					push!(totalNumExpensiveTurns[i], numExpensiveTurns[i])
					assert(numPaths[i] == length(totalPaths[i]))
					totalPaths[i][numPaths[i]], totalPaths[i][1] = totalPaths[i][1], totalPaths[i][numPaths[i]]
					totalNumExpensiveTurns[i][numPaths[i]], totalNumExpensiveTurns[i][1] = totalNumExpensiveTurns[i][1], totalNumExpensiveTurns[i][numPaths[i]]
				end
			# If we need to remove a path
			else
				worstIndex = findWorstPathIndex(totalPaths[i], totalNumExpensiveTurns[i], turnCost, newTimes)
				if worstIndex == 1
					totalPaths[i][1] = paths[i]
					totalNumExpensiveTurns[i][1] = numExpensiveTurns[i]
				else
					totalPaths[i][1], totalPaths[i][worstIndex] = paths[i], totalPaths[i][1]
					totalNumExpensiveTurns[i][1], totalNumExpensiveTurns[i][worstIndex] = numExpensiveTurns[i], totalNumExpensiveTurns[i][1]
				end
			end
		end
	end
	return totalPaths, totalNumExpensiveTurns, numPaths
end

function loadTrainingResults(directoryName::AbstractString, max_iterations::Int, manhattan::Manhattan, turnCost::Float64)
	"""
	Given the name of the directory where the algorithm output was saved, run shortest paths to determine travel times between all pairs of nodes.
	"""
	println("-- Loading results from $(directoryName)...")
	link_times = load("Outputs/$(directoryName)/manhattan-times-$(max_iterations).jld", "times")
	println("-- Computing shortest paths...")
	new_graph, new_edge_dists, new_nodes = modifyGraphForDijkstra(manhattan.network, link_times, manhattan.positions, turn_cost=turnCost)
	old_nodes = getInverseMapping(new_nodes, nv(new_graph))
	new_sp = parallelShortestPathsWithTurnsAuto(manhattan.network, new_graph, new_edge_dists, new_nodes)
	return new_sp.traveltime
end