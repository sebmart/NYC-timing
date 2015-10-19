# data_handling.jl
# New way to load data into memory to plug into Julia
# Also includes training/testing set separation and evaluation of method results

MANHATTAN_NODES = 6134

function createReducedJLDs()
	"""
	Converts reduced CSV files into JLD format with a DataFrame
	Takes in no arguments and returns nothing (designed to be run once)
	"""
	for j = 10:12
		println(j)
		f = open("../Ride_data/data_link/reduced_trip_data_$j.csv")
		df2 = DataFrame(pTime = DateTime[], dTime=DateTime[], pX=Float64[], pY=Float64[], dX=Float64[], dY=Float64[],)
		for (i, ln) in enumerate(eachline(f))
			if i == 1
				continue
			end
			if i%1000000 == 0
				println(i)
			end
			s = split(ln,",")
			p = DateTime(s[1][2:end-1], "y-m-d H:M:S")
			d = DateTime(s[2][2:end-1], "y-m-d H:M:S")
			push!(df2, [p,d, float(s[3]), float(s[4]), float(s[5]), float(strip(s[6], '\n'))])
		end
	    close(f)
	    save("../Ride_data/new_trip_data_$j.jld","df",df2)
	end
	return nothing
end

function loadRides(;startTime::Int=12, endTime::Int=14, weekdays::Bool = true, startMonth::Int=1, endMonth::Int=1, loadFromCache = true)
	"""
	Load appropriate set of rides, between startTime and endTime, startMonth and endMonth, on weekdays or weekends. Supports caching to save time.
	"""
	# Initialize output file name
	fileName = "../Ride_data/rides_$(startTime)$(endTime)_m$(startMonth)$(endMonth)"
	if weekdays
		fileName = string(fileName, "_wd.jld")
	else
		fileName = string(fileName, "_we.jld")
	end
	# Load from file if available
	if isfile(fileName)
		df2 = load(fileName, "df")
		return df2
	end
	# Create dataframe
	df2 = DataFrame()
	# Load data
	for i = startMonth:endMonth
		df = load("../Ride_data/new_trip_data_$i.jld", "df")
		if weekdays
			df2 = vcat(df2,df[((dayofweek(df[:,:pTime]) .<= 5) & (startTime .<= hour(df[:,:pTime]) .<= endTime) & (dayofweek(df[:,:dTime]) .<= 5) & (startTime .<= hour(df[:,:dTime]) .<= endTime)),:])
		else
			df2 = vcat(df2,df[((dayofweek(df[:,:pTime]) .> 5) & (startTime .<= hour(df[:,:pTime]) .<= endTime) & (dayofweek(df[:,:dTime]) .> 5) & (startTime .<= hour(df[:,:dTime]) .<= endTime)),:])
		end
	end
	# Filter out bad rides
	function goodRide(i::Int)
	    return (df2[i,:dTime] - df2[i,:pTime]).value > 60*1000 && (abs(df2[i,:pX] - df2[i,:dX]) + abs(df2[i,:pY] - df2[i,:dY])) >= 140
	end
	df2 = df2[Bool[goodRide(i) for i = 1:nrow(df2)],:]
	save(fileName, "df", df2)
	return df2
end

function splitTrainingAndTestingSets(df::DataFrame, inverseTestingSetFraction::Int)
	"""
	Takes in DataFrame with rides. Returns two DataFrames, one for training and the other for testing
	"""
	shuf = shuffle(collect(1:nrow(df)))
	testDf = df[shuf[1:round(Int,length(shuf)/inverseTestingSetFraction)],:]
	trainDf= df[shuf[(round(Int,length(shuf)/inverseTestingSetFraction) + 1):end],:]
	return trainDf, testDf
end

function loadTrainingAndTestingSets(;startTime::Int=12, endTime::Int=14, weekdays::Bool = true, startMonth::Int=1, endMonth::Int=1, loadFromCache = true)
	"""
	Loads dataframe for rides, splits off 20% for permanent testing set (NEVER to be touched), then splits training set into 80% tr-learning and 20% tr-testing sets
	"""
	fileName = "../Ride_data/rides_$(startTime)$(endTime)_m$(startMonth)$(endMonth)"
	if weekdays
		fileName = string(fileName, "_wd")
	else
		fileName = string(fileName, "_we")
	end
	testingFileName = string(fileName, "_testing.jld")
	trainLearnFileName = string(fileName, "_tr_learn.jld")
	trainTestFileName = string(fileName, "_tr_test.jld")
	if isfile(testingFileName) && isfile(trainLearnFileName) && isfile(trainTestFileName) && loadFromCache
		testDf = load(testingFileName, "df")
		trainLearnDf = load(trainLearnFileName, "df")
		trainTestDf = load(trainTestFileName, "df")
	else
		df = loadRides(startTime=startTime, endTime=endTime, weekdays=weekdays,startMonth=startMonth, endMonth=endMonth, loadFromCache=true)
		trainDf, testDf = splitTrainingAndTestingSets(df, 5)
		trainLearnDf, trainTestDf = splitTrainingAndTestingSets(trainDf, 5)
		save(testingFileName, "df", testDf)
		save(trainLearnFileName, "df", trainLearnDf)
		save(trainTestFileName, "df", trainTestDf)
	end
	return trainLearnDf, trainTestDf, testDf
end

function mapRidesToNodes(trainDf::DataFrame, nodePositions::Array{Coordinates, 1}, method::String; radius::Int = 140)
	"""
	Given the learning set of rides, an array of the node locations and the appropriate method, converts dataframe full of rides into travelTime matrix
	"""
	travelTimes = zeros(length(nodePositions), length(nodePositions))
	numRides = zeros(Int, (length(nodePositions), length(nodePositions)))
	if method == "nearest"
		# Put nodes in KDTree
		nodeTree, nodePairs = buildNodeKDTree(nodePositions)
		# Map rides to nodes, without overlap
		# Loop through rides
		for i = 1:nrow(trainDf)
			# Choose method
			if i % 50000 == 0
				println(i)
			end
			# Nearest neighbor search
			index, dist = knn(nodeTree, [trainDf[i,:pX], trainDf[i,:pY], trainDf[i,:dX], trainDf[i,:dY]], 1)
			index = index[1]
			newTime = (trainDf[i,:dTime] - trainDf[i,:pTime]).value/1000
			travelTimes[nodePairs[index,1], nodePairs[index,2]] = (travelTimes[nodePairs[index,1], nodePairs[index,2]] * numRides[nodePairs[index,1], nodePairs[index,2]] + newTime)/(numRides[nodePairs[index,1],nodePairs[index,2]]+1)
			numRides[nodePairs[index,1], nodePairs[index,2]] += 1
		end
	elseif method == "average"
		# Put rides in KDTree
		treeL = zeros(4,nrow(trainDf))
		for i in 1:4
			treeL[i,:] = trainDf[:,(i+2)]
		end
		println("--Constructing KDTree")
		@time treeL = KDTree(treeL)
		for i = 1:length(nodePositions), j=1:length(nodePositions)
			if i != j
				vec = [nodePositions[i].x, nodePositions[i].y, nodePositions[j].x, nodePositions[j].y]
				index = inball(treeL, vec, float(radius))
				for idx in index
					newTime = (trainDf[idx,:dTime] - trainDf[idx,:pTime]).value/1000
					travelTimes[i,j] = (travelTimes[i,j] * numRides[i,j] + newTime)/(numRides[i,j] + 1)
					numRides[i,j] += 1
				end
			end
		end
	end
	return travelTimes, numRides
end

function loadInputTravelTimes(nodePositions::Array{Coordinates}, method::String; startTime::Int=12, endTime::Int=14, weekdays::Bool = true, startMonth::Int=1, endMonth::Int=1, loadFromCache = true, radius::Int=140)
	"""
	Top level function, called before fast_LP to load travel times for data.
	"""
	println("**** Loading training and testing sets ****")
	# Select file name
	fileName = "../Ride_data/Input_travel_times/travel_times_$(startTime)$(endTime)_m$(startMonth)$(endMonth)"
	if weekdays
		fileName = string(fileName, "_wd")
	else
		fileName = string(fileName, "_we")
	end
	trainLearnFileName = string(fileName, "_tr_learn.jld")
	trainTestFileName = string(fileName, "_tr_test.jld")
	trainLearnDf, trainTestDf, testDf = loadTrainingAndTestingSets(startTime=startTime, endTime=endTime, weekdays=weekdays, startMonth=startMonth, endMonth=endMonth)
	if isfile(trainLearnFileName) && loadFromCache
		travelTimes = load(trainLearnFileName, "travelTimes")
		numRides = load(trainLearnFileName, "numRides")
	else
		travelTimes, numRides = mapRidesToNodes(trainLearnDf, nodePositions, method, radius=radius)
		save(trainLearnFileName, "travelTimes", travelTimes, "numRides", numRides)
	end
	if isfile(trainTestFileName) && loadFromCache
		testingTravelTimes = load(trainTestFileName, "travelTimes")
		testingNumRides = load(trainTestFileName, "numRides")
	else
		testingTravelTimes, testingNumRides = mapRidesToNodes(trainTestDf, nodePositions, method, radius=radius)
		save(trainTestFileName, "travelTimes". testingTravelTimes, "numRides", testingNumRides)
	end
	return travelTimes, numRides, trainTestDf, testingTravelTimes, testingNumRides, testDf
end

function buildNodeKDTree(nodePositions::Array{Coordinates,1})
	"""
	Given an array of node positions, returns a KDTree containing the nodePairs in 4D space. Also returns array mapping index in KDTree to pair of nodes.
	"""
	nodeTree = zeros(4, length(nodePositions)^2 - length(nodePositions))
	index = 1
	nodePairs = zeros(Int, (length(nodePositions)^2 - length(nodePositions), 2))
	for i = 1:length(nodePositions), j=1:length(nodePositions)
		if i != j
			nodeTree[1, index] = nodePositions[i].x
			nodeTree[2, index] = nodePositions[i].y
			nodeTree[3, index] = nodePositions[j].x
			nodeTree[4, index] = nodePositions[j].y
			nodePairs[index, 1] = i
			nodePairs[index, 2] = j
			index += 1
		end
	end
	@time nodeTree = KDTree(nodeTree)
	return nodeTree, nodePairs
end

function computeAverageTestingError(testingTravelTimes::Array{Float64,2}, algorithmOutput::Array{float64, 2}; num_nodes::Int = MANHATTAN_NODES)
	"""
	Given a matrix of TABs from testing Data and calculated by our algorithm, returns some error measures
	"""
	num_results = 0
	average_squared_error = 0
	average_relative_error = 0
	average_bias = 0
	for i = 1:num_nodes, j=1:num_nodes
		if testingData[i,j] > 0
			average_squared_error = (num_results * average_squared_error + (testingTravelTimes[i,j] - algorithmOutput[i,j]) ^ 2) / (num_results + 1)
			average_relative_error = (num_results * average_relative_error + abs(testingTravelTimes[i,j] - algorithmOutput[i,j])/testingTravelTimes[i,j])/(num_results + 1)
			average_bias = (num_results * average_bias + (algorithmOutput[i,j] - testingTravelTimes[i,j]))/(num_results + 1)
			num_results += 1
		end
	end
	println("------------------------------------------------")
	println("Average squared error: \t\t\t\t\t", average_squared_error)
	println("Average relative error: \t\t\t\t", average_relative_error)
	println("Average bias: \t\t\t\t\t\t", average_bias)
	println("------------------------------------------------")
	return average_squared_error, average_relative_error, average_bias

function computeTestingError(rideTestingSet::DataFrame, nodeTree::KDTree{Float64}, nodePairs::Array{Int,2}, algorithmOutput::Array{Float64,2}, method::String; radius::Int = 140)
	"""
	Uses testing set to compute error on our method.
	"""
	# Initialize vars
	num_results = 0
	average_squared_error = 0
	average_relative_error = 0
	average_bias = 0
	# Loop through testing set and choose method
	for i = 1:nrow(testDf)
		if method == "nearest"
			# Find nearest node pair to this ride
			index, dist = knn(nodeTree, [testDf[i,:pX], testDf[i,:pY], testDf[i,:dX], testDf[i,:dY]], 1)
			idx = index[1]
			# Calculate time of ride in seconds
			testingTime = (testDf[i,:dTime] - testDf[i,:pTime]).value/1000 # don't forget to convert
			# Calculate predicted time
			predictedTime = algorithmOutput[nodePairs[idx,1],nodePairs[idx,2]]
			# Update errors
			if predictedTime != 0
				average_squared_error = (num_results * average_squared_error + (testingTime - predictedTime) ^ 2) / (num_results + 1)
				average_relative_error = (num_results * average_relative_error + abs(testingTime - predictedTime)/testingTime) / (num_results + 1)
				average_bias = (num_results * average_bias + (predictedTime - testingTime))/(num_results + 1)
				num_results += 1
			end
		elseif method == "average"
			index = inball(nodeTree, [testDf[i,:pX], testDf[i,:pY], testDf[i,:dX], testDf[i,:dY]], float(radius))
			# Compute testing time in seconds
			testingTime = (testDf[i,:dTime] - testDf[i,:pTime]).value/1000
			# Compute predicted time by averaging neighboring node pairs (within radius)
			predictedTime = 0.
			numPoints = 0
			for idx in index
				predictedTime = (predictedTime * numPoints + algorithmOutput[nodePairs[idx,1],nodePairs[idx,2]])/(numPoints + 1)
				numPoints += 1
			end
			# Update errors (if any neighbors were found within 140 meters)
			if predictedTime != 0
				average_squared_error = (num_results * average_squared_error + (testingTime - predictedTime) ^ 2) / (num_results + 1)
				average_relative_error = (num_results * average_relative_error + abs(testingTime - predictedTime)/testingTime) / (num_results + 1)
				average_bias = (num_results * average_bias + (predictedTime - testingTime))/(num_results + 1)
				num_results += 1
			end
		end
	end
	println("------------------------------------------------")
	println("Average squared error: \t\t\t\t\t", average_squared_error)
	println("Average relative error: \t\t\t\t", average_relative_error)
	println("Average bias: \t\t\t\t\t\t", average_bias)
	println("------------------------------------------------")
	return average_squared_error, average_relative_error, average_bias
end