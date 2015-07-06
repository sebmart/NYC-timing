cd("/Users/bzeng/Dropbox (MIT)/7 Coding/UROP/NYC-timing")
rides = readcsv("manhattan_rides_wd_1214.csv")
numRides = convert(Int, length(rides) / 5)
originalIndices = shuffle(collect(1:numRides))
splitIndex = convert(Int, floor(length(originalIndices) * 8/10))
trainingIndices = shuffle(originalIndices[1:splitIndex])
testingIndices = originalIndices[splitIndex + 1:length(originalIndices)]

k = 5
trainingSubsetSize = convert(Int, floor(size(trainingIndices)[1] / k)) 
trainingSubsets = reshape(trainingIndices[1: k * trainingSubsetSize], k, trainingSubsetSize)

# Based on the value of k, partitions the trainingSet into further trainingSubsets
function partitionSet(k::Int64)
	subsets = collect(1:k)
	partitions = Int64[]
	for i = 1:k
		if i == 1
			partitions = vcat(partitions, subsets[2:k])
		elseif i == k
			partitions = vcat(partitions, subsets[1:k - 1])
		else
			partitions = vcat(partitions, vcat(subsets[1:i - 1], subsets[i + 1:k]))
		end
	end
	return reshape(partitions, k, k-1)
end

# Given a set of parameters and a set of rides, determines hat(T_AB) (estimated road times)
function Arthur(paramaters::Float64, indices::Vector{Int64})
	return 0
end

# Computes the test error between expected road times and computed road times
function testError(computedRoadTime::Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}, originalRoadTime::Base.SparseMatrix.SparseMatrixCSC{Float64,Int64})
	dimension = convert(Int, length(computedRoadTime) ^ 0.5)
	bias = 0
	variance = 0
	for i = 1:dimension
		for j = 1:dimension
			difference = computedRoadTime[i, j] - originalRoadTime[i, j]
			bias += difference
			variance += difference ^ 2
		end
	end
	return (bias, variance)
end

partitions = partitionSet(k)
bias = 0
variance = 0
for i = 1:k
	subsets = partitions[i, 1:k-1]
	test = Int64[]
	for subset in subsets
		for j = 1:trainingSubsetSize
			push!(test, trainingSubsets[subset, j])
		end
	end
	computedRoadTime = Arthur(paramater, test)
	error = testError(computedRoadTime, originalRoadTime)
	bias += error[1]
	variance += error[2]
end

