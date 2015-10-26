# runTest.jl
# Compares matrix of TABs generated from the algorithm to real testing data
# Authored by Arthur J Delarue on 7/8/15

MANHATTAN_NODES = 6134

function compute_error(testingData::AbstractArray{Float64, 2}, trainingData::AbstractArray{Float64, 2}; num_nodes::Int = MANHATTAN_NODES)
	"""
	Given a matrix of TABs from testing Data and calculated by our algorithm, returns some error measures
	"""
	num_results = 0
	average_squared_error = 0
	average_relative_error = 0
	average_bias = 0
	for i = 1:num_nodes, j=1:num_nodes
		if testingData[i,j] > 0
			average_squared_error = (num_results * average_squared_error + (testingData[i,j] - trainingData[i,j]) ^ 2) / (num_results + 1)
			average_relative_error = (num_results * average_relative_error + abs(testingData[i,j] - trainingData[i,j])/testingData[i,j])/(num_results + 1)
			average_bias = (num_results * average_bias + (trainingData[i,j] - testingData[i,j]))/(num_results + 1)
			num_results += 1
		end
	end
	println("------------------------------------------------")
	println("Average squared error: \t\t\t\t\t", average_squared_error)
	println("Average relative error: \t\t\t\t", average_relative_error)
	println("Average bias: \t\t\t\t\t\t", average_bias)
	println("------------------------------------------------")
	return average_squared_error, average_relative_error, average_bias
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


# trainingData = loadTrainingResults(NAME, 5, manhattan, 10.0)
# compute_error(testingData, trainingData)