using HDF5, JLD, KDTrees, DataFrames, Base.Dates

function selectingTimes (start::Int64, finish::Int64, radius::Float64)
	cd("/Users/bzeng/Dropbox (MIT)/7 Coding/UROP/Data")
	start = 12
	finish = 14
	startingHour = start
	endingHour = finish
	radius = 140.0
	nodePositions = load("manhattan.jld","positions")
	numRides = zeros(Int64, length(nodePositions), length(nodePositions))
	avgTime = zeros(Float32, length(nodePositions), length(nodePositions))
	
	nodes = Array(Float64, (2,length(nodePositions)))
	for (i,p) in enumerate(nodePositions)
	    nodes[1,i] = p[1]
	    nodes[2,i] = p[2]
	end

	tree = KDTree(nodes)

	function coordinatesToNode(x::Float64, y::Float64)
		return inball(tree, [x, y], radius, false)
	end

	for i = 1:12
		println(i)
		df = 0
		@time df = readtable("reduced_trip_data_$(i).csv");
		for j = 1:nrow(df)
			if j % 1000 == 0
				print("\r $(j * 100 / nrow(df))                   ")
			end
			pickup = df[j, :pickup_datetime]
			dropoff = df[j, :dropoff_datetime]
			pickupTime = DateTime(pickup, "y-m-d H:M:S")
			dropoffTime = DateTime(dropoff, "y-m-d H:M:S")
			pickupFlag = 1 <= dayofweek(pickupTime) <= 5 && startingHour <= hour(pickupTime) <= endingHour
			dropoffFlag = 1 <= dayofweek(dropoffTime) <= 5 && startingHour <= hour(dropoffTime) <= endingHour
			if pickupFlag && dropoffFlag
				nodeSrcList = coordinatesToNode(df[j, :pickup_x], df[j, :pickup_y])
				nodeDstList = coordinatesToNode(df[j, :dropoff_x], df[j, :dropoff_y])
				for nodeSrc in nodeSrcList, nodeDst in nodeDstList
					numRides[nodeSrc, nodeDst] += 1
					prevTime = avgTime[nodeSrc, nodeDst]
					newTime = (dropoffTime - pickupTime).value / 1000.0 
					avgTime[nodeSrc, nodeDst] = (newTime + prevTime * (numRides[nodeSrc, nodeDst] - 1)) / numRides[nodeSrc, nodeDst]
				end
			end
		end
	end

	df = DataFrame(node1 = Int64[], node2 = Int64[], numberOfRides = Int64[], averageTime = Float32[])
	for src = 1:length(nodePositions), dst = 1:length(nodePositions)	
		if numRides[src, dst] > 0
			push!(df, [src dst numRides[src, dst] avgTime[src, dst]])
		end
	end

	@time writetable("trip_data_$(startingHour)_$(endingHour).csv", df)
end

@time selectingTimes(12, 14, 140.0)

