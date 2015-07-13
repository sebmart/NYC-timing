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
	return average_squared_error, average_relative_error, average_bias
end

function loadTrainingResults(directoryName::String, max_iterations::Int, manhattan::Manhattan, turnCost::Float64)
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

# manhattan = loadCityGraph()
# NAME = "new_r140_minr3_i5_wd_1214_avg_clust50_rides50000"
# trainingData = loadTrainingResults(NAME, 2, manhattan, 10.0)
# testingData, numRides = loadNewTravelTimeData(trainOrTest="testing", radius = 140, times = "1214", preprocess = false)
# compute_error(testingData, trainingData)