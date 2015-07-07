using HDF5, JLD, KDTrees, DataFrames, Base.Dates
cd("/Users/bzeng/Dropbox (MIT)/7 Coding/UROP/Data")
function processTestingData(radius::Float64, starting::Int64, ending::Int64)
	testingDF = readtable("testing_r$(radius)_wd_$(starting)$(ending).csv")
	nodePositions = load("manhattan.jld","positions")
	roadTimes = zeros(Int64, length(nodePositions), length(nodePositions))
	for i = 1:nrow(testingDF)
		roadTimes[testingDF[i, :node1], testingDF[i, :node2]] = testingDF[i, :averageTime]
	end
	return roadTimes
end
