# averaging.jl
# Small module to test lower bounds on error in Manhattan
# Authored by Arthur Delarue on 11/30/15

function average_manhattan_results(fileList)
	times = spzeros(6134, 6134)
	for (i,fileName) in enumerate(fileList)
		if i == 1
			times = load(fileName, "times")
		else
			tmp_times = load(fileName, "times")
			times = (times .* (i - 1) + tmp_times) ./ i
		end
	end
	return times
end

manhattan = loadCityGraph()
graph = manhattan.network
travel_times, num_rides, trainTestDf, testing_travel_times, testing_num_rides, testDf = loadInputTravelTimes(manhattan.positions, "average", year=2013, startTime=12, endTime=14, startMonth=1, endMonth=12, radius = 200)
nodeTree, nodePairs = buildNodeKDTree(manhattan.positions)

newTimes = average_manhattan_results([
	"Outputs/fast_r200_minr1_i50_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st50000_add0_thr0.7_rem0_thr0.3_every1_gcu_tc0.0/manhattan-times-31.jld",
	"Outputs/fast_r200_minr1_i50_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st50000_add0_thr0.7_rem0_thr0.3_every2_gcu_tc0.0/manhattan-times-50.jld", 
	"Outputs/fast_r200_minr1_i50_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st50000_add0_thr0.7_rem0_thr0.3_every3_gcu_tc0.0/manhattan-times-50.jld",
	"Outputs/fast_r200_minr1_i50_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st50000_add0_thr0.7_rem0_thr0.3_every4_gcu_tc0.0/manhattan-times-50.jld"])

new_graph, new_edge_dists, new_nodes, new_edge_isExpensive = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=0.0)
old_nodes = getInverseMapping(new_nodes, nv(new_graph))
new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
avg_sq_err, avg_rel_err, avg_bias = computeTestingError(trainTestDf, nodeTree, nodePairs, new_sp.traveltime, "average")