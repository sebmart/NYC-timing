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
	numNodePairsToRemove::Int = 0)
	newSrcs = deepcopy(srcs)
	newDsts = deepcopy(dsts)
	newPairs = deepcopy(pairs)
	newTotalPaths = deepcopy(totalPaths)
	newTotalNumExpensiveTurns = deepcopy(totalNumExpensiveTurns)
	newNumPaths = deepcopy(numPaths)

	# Create lookup set for constraints already under consideration
	nodePairsAlreadyIn = Set([(srcs[i], dsts[i]) for i = 1:length(srcs)])
	# Remove worst constraints
	# TODO
	# Select new constraints
	indicesVector = Tuple{Int,Int}[]
	sizehint!(indicesVector, 12000000)
	for i = 1:size(travelTimes)[1], j = 1:size(travelTimes)[1]
		if travelTimes[i,j] > 0
			push!(indicesVector, (i,j))
		end
	end
	# compute distance between calculated times and data
	errorVector = zeros(length(indicesVector))
	for i = 1:length(indicesVector)
		errorVector[i] = sqrt(numRides[indicesVector[i][1], indicesVector[i][2]]/travelTimes[indicesVector[i][1], indicesVector[i][2]]) * abs(travelTimes[indicesVector[i][1], indicesVector[i][2]] - calculatedTravelTimes[indicesVector[i][1], indicesVector[i][2]])
	end
	# sort
	p = sortperm(errorVector)
	startIndex = int(round(0.9 * length(indicesVector)))
	i = 0
	count = 0
	while i < startIndex - 1 && count < numNodePairsToAdd
		if !(indicesVector[p[startIndex - i]] in nodePairsAlreadyIn)
			push!(newSrcs, indicesVector[p[startIndex - i]][1])
			push!(newDsts, indicesVector[p[startIndex - i]][2])
			push!(newTotalPaths, Array{Int}[])
			push!(newTotalNumExpensiveTurns, Int[])
			push!(newNumPaths, 0)
			push!(newPairs[indicesVector[p[startIndex - i]][1]], indicesVector[p[startIndex - i]][2])
			sort!(newPairs[indicesVector[p[startIndex - i]][1]])
			count += 1
		end
		i += 1
	end
	return newSrcs, newDsts, newTotalPaths, newTotalNumExpensiveTurns, newNumPaths, newPairs
end

function generate_indices(totalEntries::Int, selectedEntries::Int, percentile::Float64)
	middle = percentile * totalEntries
	low = int(round(middle - selectedEntries/2))
	high = low + selectedEntries - 1
	return low:high
end