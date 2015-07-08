using HDF5, JLD, KDTrees, DataFrames
using Dates

function testingTimes (start::Int64, finish::Int64, radius::Float64)
	startingHour = start
	endingHour = finish
	radius = 140.0
	nodePositions = load("../Cities/Manhattan/manhattan.jld","positions")

	trainingNumRides = zeros(Int64, length(nodePositions), length(nodePositions))
	trainingAvgTime = zeros(Float32, length(nodePositions), length(nodePositions))
	testingNumRides = zeros(Int64, length(nodePositions), length(nodePositions))
	testingAvgTime = zeros(Float32, length(nodePositions), length(nodePositions))

	nodes = Array(Float64, (2,length(nodePositions)))
	for (i,p) in enumerate(nodePositions)
	    nodes[1,i] = p[1]
	    nodes[2,i] = p[2]
	end

	tree = KDTree(nodes)

	function coordinatesToNode(x::Float64, y::Float64)
		return inball(tree, [x, y], radius, false)
	end

	function computeDistance(px::Float64, py::Float64, dx::Float64, dy::Float64)
		return sqrt((px - dx) ^ 2 + (py - dy) ^ 2)
	end

	for i = 1:12
		println(i)
		df = 0
		@time df = readtable("../../Data/reduced_trip_data_$(i).csv");
		for j = 1:nrow(df)
			if j % 10000 == 0
				print("\r $(j * 100 / nrow(df))%            ")
			end
			pickup = df[j, :pickup_datetime]
			dropoff = df[j, :dropoff_datetime]
			pickupTime = DateTime(pickup, "y-m-d H:M:S")
			dropoffTime = DateTime(dropoff, "y-m-d H:M:S")
			pickupFlag = 1 <= dayofweek(pickupTime) <= 5 && startingHour <= hour(pickupTime) <= endingHour
			dropoffFlag = 1 <= dayofweek(dropoffTime) <= 5 && startingHour <= hour(dropoffTime) <= endingHour
			if pickupFlag && dropoffFlag
				if computeDistance(df[j, :pickup_x], df[j, :pickup_y], df[j, :dropoff_x], df[j, :dropoff_y]) > 140.0 
					nodeSrcList = coordinatesToNode(df[j, :pickup_x], df[j, :pickup_y])
					nodeDstList = coordinatesToNode(df[j, :dropoff_x], df[j, :dropoff_y])
					for nodeSrc in nodeSrcList, nodeDst in nodeDstList
						if rand() > 0.2
							trainingNumRides[nodeSrc, nodeDst] += 1
							prevTime = trainingAvgTime[nodeSrc, nodeDst]
							newTime = (dropoffTime - pickupTime).value / 1000.0
							trainingAvgTime[nodeSrc, nodeDst] = (newTime + prevTime * (trainingNumRides[nodeSrc, nodeDst] - 1)) / trainingNumRides[nodeSrc, nodeDst]
						else
							testingNumRides[nodeSrc, nodeDst] += 1
							prevTime = testingAvgTime[nodeSrc, nodeDst]
							newTime = (dropoffTime - pickupTime).value / 1000.0 
							testingAvgTime[nodeSrc, nodeDst] = (newTime + prevTime * (testingNumRides[nodeSrc, nodeDst] - 1)) / testingNumRides[nodeSrc, nodeDst]
						end
					end
				end
			end
		end
	end

	trainingDF = DataFrame(node1 = Int64[], node2 = Int64[], numberOfRides = Int64[], averageTime = Float32[])
	testingDF = DataFrame(node1 = Int64[], node2 = Int64[], numberOfRides = Int64[], averageTime = Float32[])
	for src = 1:length(nodePositions), dst = 1:length(nodePositions)	
		if trainingNumRides[src, dst] > 0
			push!(trainingDF, [src dst trainingNumRides[src, dst] trainingAvgTime[src, dst]])
		end
		if testingNumRides[src, dst] > 0
			push!(testingDF, [src dst testingNumRides[src, dst] testingAvgTime[src, dst]])
		end
	end
	
	startingString = string(startingHour)
	endingString = string(endingHour)
	if startingHour < 10
		startingString = "0" * startingString
	end
	if endingHour < 10
		endingString = "0" * endingString
	end
	writetable("../../Travel_times/training_r$(int(radius))_wd_$(startingString)$(endingString).csv", trainingDF)
	writetable("../../Travel_times/testing_r$(int(radius))_wd_$(startingString)$(endingString).csv", testingDF)
end

# testingTimes(12, 14, 140.0)
# testingTimes(2, 4, 140.0)
# testingTimes(7, 9, 140.0)
