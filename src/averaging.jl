# averaging.jl
# Small module to test lower bounds on error in Manhattan
# Authored by Arthur Delarue on 11/30/15

manhattan = loadCityGraph()
graph = manhattan.network
travel_times, num_rides, trainTestDf, testing_travel_times, testing_num_rides, testDf = loadInputTravelTimes(manhattan.positions, "average", year=2013, startTime=12, endTime=14, startMonth=1, endMonth=12, radius = 200)
nodeTree, nodePairs = buildNodeKDTree(manhattan.positions)

function average_manhattan_results(fileList)
	travelTimes = zeros(6134, 6134)
	for (i,fileName) in enumerate(fileList)
		newTimes = load(fileName, "times")
		new_graph, new_edge_dists, new_nodes, new_edge_isExpensive = modifyGraphForDijkstra(graph, newTimes, manhattan.positions, turn_cost=0.0)
		old_nodes = getInverseMapping(new_nodes, nv(new_graph))
		new_sp = parallelShortestPathsWithTurnsAuto(graph, new_graph, new_edge_dists, new_nodes)
		if i == 1
			travelTimes = new_sp.traveltime
		else
			travelTimes = (travelTimes .* (i - 1) + new_sp.traveltime) ./ i
		end
	end
	return travelTimes
end

newTravelTimes = average_manhattan_results([
	"Outputs/fast_r200_minr1_i50_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st50000_add0_thr0.7_rem0_thr0.3_every1_gcu_tc0.0/manhattan-times-31.jld",
	"Outputs/fast_r200_minr1_i50_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st50000_add0_thr0.7_rem0_thr0.3_every2_gcu_tc0.0/manhattan-times-50.jld", 
	"Outputs/fast_r200_minr1_i50_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st50000_add0_thr0.7_rem0_thr0.3_every3_gcu_tc0.0/manhattan-times-50.jld",
	"Outputs/fast_r200_minr1_i50_wd_1214_m112_y2013_relaxed_smooth_ppc3_dataaverage_errorboth_dynConstr_st50000_add0_thr0.7_rem0_thr0.3_every4_gcu_tc0.0/manhattan-times-50.jld"])


avg_sq_err, avg_rel_err, avg_bias = computeTestingError(trainTestDf, nodeTree, nodePairs, newTravelTimes, "average")