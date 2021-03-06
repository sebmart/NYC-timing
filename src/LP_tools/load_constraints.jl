# load_constraints.jl
# Gets city information from manhattan.jld object and travel time information from Travel time CSV
# Authored by Arthur J Delarue on 6/8/15

function saveManhattan(pb::TaxiProblem, name::AbstractString; compress=false)
  save("Cities/Saved/$name.jld", "pb", pb, compress=compress)
end

function loadManhattan(name::AbstractString)
  pb = load("Cities/Saved/$name.jld","pb")
  return pb
end

function saveRoadTimes(times::SparseMatrixCSC{Float64,Int}, name::AbstractString; compress=false)
  save("Outputs/$name.jld", "times", times, compress=compress)
end

function loadRoadTimes(name::AbstractString)
  times = load("Outputs/$name.jld","times")
  return times
end

MANHATTAN_NODES = 6134
# Speed limits in km/h from OpenStreetMap/src/speeds.jl
URBAN_SPEEDS = [95,72,48,32,22,12,8,5]

function loadCityGraph(;fromScratch::Bool=false, useShortestPaths::Bool=true,fixSpeeds=false)
	"""
	Load Manhattan object and return it
	"""
	println("**** Loading city graph ****")
	if fromScratch
		manhattan = Manhattan(sp=useShortestPaths)
		# if useShortestPaths
		# 	sp = manhattan.sp
		# else
		# 	sp = parallelShortestPathsAuto(manhattan.network, manhattan.roadTime, manhattan.roadCost)
		# 	manhattan.sp = sp
		# end
		if fixSpeeds
			# Make sure the speeds are in URBAN_SPEEDS
			# If not match speed to closest actual speed
			for i=1:MANHATTAN_NODES,j=1:MANHATTAN_NODES
				if manhattan.distances[i,j] != 0
					# Don't forget to convert from m/s to km/h
					if manhattan.roadTime[i,j] != 0
						speed = 3.6 * manhattan.distances[i,j]/manhattan.roadTime[i,j]
					else
						speed = 3.6 * manhattan.distances[i,j]/0.1
					end
					road_type = indmin(abs(URBAN_SPEEDS - speed))
					if manhattan.roadTime[i,j] == 0
						println(manhattan.roadTime[i,j], "\t", manhattan.distances[i,j], "\t", speed, "\t", URBAN_SPEEDS[road_type])
					end
					manhattan.roadTime[i,j] = manhattan.distances[i,j]/URBAN_SPEEDS[road_type]
				end
			end
			shortestPaths!(manhattan)
		end
		saveManhattan(manhattan, "full_manhattan")
	else
		manhattan = loadManhattan("full_manhattan")
	end
	return manhattan
end

function loadTravelTimeData(;radius::Int=40, times::AbstractString= "1214", min_rides::Int=4, num_nodes::Int=MANHATTAN_NODES, preprocess::Bool=false, num_clusters::Int64 = 50, sampleSize::Int64=10000) 
	# Hardcoded for simplicity. Map of Manhattan not likely to change anytime soon
	"""
	Return array of travel times between pairs of nodes A,B obtained from NYC dataset
	"""
	println("**** Loading travel times ****")
	# Load Manhattan highway nodes
	highwayNodes = load("Cities/Saved/highwayNodes.jld", "highwayNodes")
	if preprocess
		dataFile = "../Travel_times/travel_times_r$(radius)_wd_$(times)_clust$(num_clusters)_rides$(sampleSize)_minr$(min_rides).csv"
	else
		dataFile = "../Travel_times/travel_times_r$(radius)_wd_$(times).csv"
	end
	println("-- Loading from $dataFile...")
	data = readcsv(dataFile)
	nodePairs = data[:,1]
	average = data[:,2]
	minRides = data[:,5]

	# Put travel time data in an array
	travelTimes = zeros(num_nodes,num_nodes)
	for (i, element) in enumerate(nodePairs)
		node = split(element, ";")
		# Remember +1 because Python indexes from 0 but Julia from 1
		src = int(node[1]) + 1
		dest = int(node[2]) + 1
		if minRides[i] >= min_rides && average[i] >= 60 && src != dest && !(src in highwayNodes) && !(dest in highwayNodes)
			travelTimes[src, dest] = average[i]
		end
	end
	return travelTimes
end

function loadNewTravelTimeData(;trainOrTest::AbstractString="training", radius::Int=40, times::AbstractString="1214", num_nodes::Int=MANHATTAN_NODES, preprocess::Bool=false, num_clusters::Int64=50, sampleSize::Int64=10000, min_rides::Int=4, saveTestingMatrix::Bool=false, loadTestingMatrixDirectly=false)
	"""
	Return array of travel times between pairs of nodes A,B obtained from NYC dataset, as well as number of rides for each node pair.
	"""
	println("**** Loading $(trainOrTest) data ****")
	# Load Manhattan highway nodes
	highwayNodes = load("Cities/Saved/highwayNodes.jld", "highwayNodes")
	if preprocess
		dataFile = "../Travel_times/$(trainOrTest)_r$(radius)_wd_$(times)_clust$(num_clusters)_rides$(sampleSize)_minr$(min_rides).csv"
	else
		dataFile = "../Travel_times/$(trainOrTest)_r$(radius)_wd_$(times).csv"
	end
	# Special case: load testing matrix directly
	if preprocess == false && loadTestingMatrixDirectly == true
		println("-- Loading directly from JLD archive")
		travelTimes = load("../Travel_times/$(trainOrTest)_r$(radius)_wd_$(times).jld", "travelTimes")
		numRides = load("../Travel_times/$(trainOrTest)_r$(radius)_wd_$(times).jld", "numRides")
		return travelTimes, numRides
	end
	println("-- Loading from $dataFile...")
	data = readcsv(dataFile)
	nodePairs = data[:,1]
	averageTime = data[:,2]
	minRides = data[:,5]
	travelTimes = zeros(num_nodes,num_nodes)
	numRides = zeros(Int, (num_nodes, num_nodes))
	for (i, element) in enumerate(nodePairs)
		node = split(element, ";")
		src = int(node[1]) + 1
		dest = int(node[2]) + 1
		if src != dest && !(src in highwayNodes) && !(dest in highwayNodes) && averageTime[i] >= 60
			travelTimes[src, dest] = averageTime[i]
			numRides[src, dest] = minRides[i]
		end
	end
	if preprocess == false && saveTestingMatrix == true
		println("-- Saving to JLD file for speed")
		save("../Travel_times/$(trainOrTest)_r$(radius)_wd_$(times).jld", "travelTimes", travelTimes, "numRides", numRides)
	end
	return travelTimes, numRides
end

function chooseConstraints(travelTimes::Array{Float64, 2}, numRides::Array{Int, 2}; num_nodes::Int = MANHATTAN_NODES, sample_size::Int = 5000)
	"""
	Given a matrix of travel times and number of rides, randomly selects sample_size nonzero entries
	"""
	newTravelTimes = zeros(num_nodes, num_nodes)
	newNumRides = zeros(Int, (num_nodes, num_nodes))
	for i = 1:sample_size
		startNode = rand(1:num_nodes)
		endNode = rand(1:num_nodes)
		while true
			# Keep them if not added yet
			if numRides[startNode, endNode] > 0 && newNumRides[startNode, endNode] == 0
				break
			end
			startNode = rand(1:num_nodes)
			endNode = rand(1:num_nodes)
		end
		newNumRides[startNode, endNode] = numRides[startNode, endNode]
		newTravelTimes[startNode, endNode] = travelTimes[startNode, endNode]
	end
	return newTravelTimes, newNumRides
end