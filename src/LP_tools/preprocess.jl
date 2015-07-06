# preprocess.jl
# Runs clustering algorithms to harmonize set of constraints over New York City
# Authored on 6/22/15 by Brandon Zeng

function outputPreprocessedConstraints(manhattan::Manhattan;radius::Int=40, numClusters::Int64=50, minRides::Int64=4, sampleSize::Int64=5000, overwrite::Bool=false)
	"""
	Checks if particular preprocessed file already exists. If so, does nothing. Otherwise, generates it.
	"""
	if overwrite || !isfile("../Travel_times/travel_times_r$(radius)_wd_1214_clust$(numClusters)_rides$(sampleSize)_minr$(minRides).csv")
		println("**** Selecting travel times ****")
		inputName = "../Travel_times/travel_times_r$(radius)_wd_1214"
		preprocess(inputName, numClusters, minRides, sampleSize, manhattan)
	else
		println("**** Selected travel times found ****")
	end
end

function preprocess(inputName::String, numClusters::Int64, minRides::Int64, sampleSize::Int64, manhattan::Manhattan)
	"""
	Clusters the nodes in Manhattan into numClusters clusters. Then reduces the number of constraints to sampleSize, so that the constraints cover as many cluster pairs as possible.
	"""
	println("---- Clustering ----")

	# Load list of nodes that are part of a highway (to exclude rides with these nodes)
	highwayNodes = load("Cities/Saved/highwayNodes.jld", "highwayNodes")

	# Load Manhattan and cluster nodes geographically
	man = manhattan
	coordinatesArray = Float64[]
	for c in man.positions
		push!(coordinatesArray, c.x)
		push!(coordinatesArray, c.y)
	end
	coordinateMatrix = reshape(coordinatesArray, 2, 6134)
	R = kmeans(coordinateMatrix, numClusters, maxiter = 1000)

	# Read in rides
	println("---- Reading CSV ----")
	rides = readcsv("$(inputName).csv")
	nodes = rides[:,1]
	average = rides[:,2]
	stddev = rides[:,3]
	stderr = rides[:,4]
	numRides = rides[:,5]

	function extractNodes(nodePair::SubString{ASCIIString})
		source = parse(Int, split(nodePair, ";")[1]) + 1
		dest = parse(Int, split(nodePair, ";")[2]) + 1
		return source, dest
	end

	function classify(source::Int, dest::Int, clusterAssignments::Vector{Int64})
		return (clusterAssignments[source], clusterAssignments[dest])
	end

	givenRides = Dict{Tuple{Int64,Int64}, Vector{Int}}()
	selectedRides = Dict{Tuple{Int64,Int64}, Vector{Bool}}()

	println("---- Loading rides ----")
	# Load rides into memory
	for i = 1:length(nodes)
		if i % 10000 == 0
			@printf("                                \r")
			@printf("Progress: %.1f%% completed\r", 100*i/length(nodes))
		end
		# Check that there are enough rides and that the time is not ridiculously small
		if numRides[i] >= minRides && average[i] >= 60
			src, dest = extractNodes(nodes[i])
			# Check that we are not starting or ending on a highway (makes no sense)
			if !(src in highwayNodes) && !(dest in highwayNodes) && src != dest
				clusterPair = classify(src, dest, R.assignments)
				if haskey(givenRides, clusterPair)
					push!(givenRides[clusterPair], i)
					push!(selectedRides[clusterPair], false)
				else
					givenRides[clusterPair] = [i]
					selectedRides[clusterPair] = [false]
				end
			end
		end
	end

	# Print out some info about how many rides we have
	lengths = [length(givenRides[key]) for key in keys(givenRides)]
	println("Number of cluster pairs with data: ", length(givenRides))
	println("Average number of rides/cluster pair: ", length(nodes)/length(givenRides))
	println("Median number of rides/cluster pair: ", median(lengths))
	totalSelections = sampleSize
	clusterSelection = round(Int, totalSelections / length(givenRides)) + 1
	result = Int[]

	println("---- Selecting rides ----")
	# Select number of rides from each cluster ("cut out the fat")
	for key in keys(givenRides)
		if length(givenRides[key]) <= clusterSelection
			result = vcat(result, givenRides[key])
		else
			count = 0
			while count < clusterSelection
				index = rand(1:length(givenRides[key]))
				if !selectedRides[key][index]
					selectedRides[key][index] = true
					push!(result, givenRides[key][index])
					count += 1
				end
			end
		end
	end

	# Sort the ride ids thus obtained so that we can quickly extract them
	sort!(result)
	writecsv("$(inputName)_clust$(numClusters)_rides$(sampleSize)_minr$(minRides).csv", rides[result,:])
end