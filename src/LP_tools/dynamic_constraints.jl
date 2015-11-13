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

	newSrcs = deepcopy(srcs)
	newDsts = deepcopy(dsts)
	newPairs = deepcopy(pairs)
	newTotalPaths = deepcopy(totalPaths)
	newTotalNumExpensiveTurns = deepcopy(totalNumExpensiveTurns)
	newNumPaths = deepcopy(numPaths)

	# Create lookup set for constraints already under consideration
	nodePairsAlreadyIn = Set([(srcs[i], dsts[i]) for i = 1:length(srcs)])
	# Select new constraints
	indicesVector = Tuple{Int,Int}[]
	sizehint!(indicesVector, 12000000)
	for i = 1:size(travelTimes)[1], j = 1:size(travelTimes)[1]
		if travelTimes[i,j] > 0
			push!(indicesVector, (i,j))
		end
	end
	# compute distance between calculated times and data
	newErrorVector = zeros(length(indicesVector))
	for i = 1:length(indicesVector)
		newErrorVector[i] = sqrt(numRides[indicesVector[i][1], indicesVector[i][2]]/travelTimes[indicesVector[i][1], indicesVector[i][2]]) * abs(travelTimes[indicesVector[i][1], indicesVector[i][2]] - calculatedTravelTimes[indicesVector[i][1], indicesVector[i][2]])
	end
	errorVector = zeros(length(srcs))
	for i = 1:length(srcs)
		errorVector[i] = sqrt(numRides[srcs[i], dsts[i]]/travelTimes[srcs[i], dsts[i]]) * abs(travelTimes[srcs[i], dsts[i]] - calculatedTravelTimes[srcs[i], dsts[i]])
	end
	# sort
	p2 = sortperm(errorVector)
	p = sortperm(newErrorVector)
	minIndexToAdd = round(Int, addOnlyIfAbovePercentile * length(indicesVector))
	maxIndexToRemove = round(Int, removeOnlyIfBelowPercentile * length(newSrcs))
	# Remove best node pairs
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
	newSrcs = newSrcs[indicesToKeep]
	newDsts = newDsts[indicesToKeep]
	newTotalNumExpensiveTurns = newTotalNumExpensiveTurns[indicesToKeep]
	newNumPaths = newNumPaths[indicesToKeep]
	newTotalPaths = newTotalPaths[indicesToKeep]
	# Add poor node pairs
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
	paths::Array{Any, 1},
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
	if dynamicConstraints && globalConstraintUpdate
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
		if nPaths < maxNumPathsPerOD * length(totalPaths)
			errorVector = zeros(nPaths)
			indicesVector = Tuple{Int,Int,Int,Int}[]
			sizehint!(indicesVector, nPaths)
			for i = 1:length(totalPaths), j=1:numPaths[i]
				if j > 1
					push!(indicesVector, (i,j,srcs[i],dsts[i]))
				end
			end
			for i=1:length(indicesVector)
				errorVector[i] = abs(sum([times[totalPaths[indicesVector[i][1]][indicesVector[i][2]][a], totalPaths[indicesVector[i][1]][indicesVector[i][2]][a+1]] for a = 1:(length(totalPaths[indicesVector[i][1]][indicesVector[i][2]])-1)]) + totalNumExpensiveTurns[indicesVector[i][1]][indicesVector[i][2]] * turnCost - sum([times[totalPaths[indicesVector[i][1]][1][a], totalPaths[indicesVector[i][1]][1][a+1]] for a = 1:(length(totalPaths[indicesVector[i][1]][1])-1)]) - totalNumExpensiveTurns[indicesVector[i][1]][1] * turnCost)/(travelTimes[indicesVector[i][3], indicesVector[i][4]])
			end
			p = sortperm(errorVector)
			numPathsToRemove = maxNumPathsPerOD * length(totalPaths) - nPaths
			pathsToKeep = [collect(1:numPaths[i]) for i=1:length(numPaths)]
			for i = 1:numPathsToRemove
				idx = findfirst(pathsToKeep[indicesVector[p[end+1-i]][1]], indicesVector[p[end+1-i]][2])
				splice!(pathsToKeep[indicesVector[p[end+1-i]][1]], idx)
			end
			totalPaths = [[totalPaths[i][j] for j in pathsToKeep[i]] for i = 1:length(totalPaths)]
			totalNumExpensiveTurns = [[totalNumExpensiveTurns[i][j] for j in pathsToKeep[i]] for i = 1:length(totalNumExpensiveTurns)]
			numPaths = [length(pathsToKeep[i]) for i = 1:length(numPaths)]
		end
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