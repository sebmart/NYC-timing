# algorithmOutput.jl
# Use manhattanVisualizations to see step-by-step progression of algorithm
# Authored by Arthur J Delarue on 6/16/15

using SFML, LightGraphs, JLD
using Base.Collections

include("../src/LP_tools/definitions.jl")
include("../src/LP_tools/manhattan.jl")

include("manhattanVisualization.jl")

DIRECTORY = "new_r140_minr1_i5_wd_1214_avgtr_clust50_rides50000"
MAX_ROUNDS = 5

function showAlgorithmEvolution(dir=DIRECTORY, firstRound::Int=1, secondRound::Int=2)
	man = Manhattan()
	oldTimes = load("/Users/arthurdelarue/Desktop/VRP-remote/Outputs/$(dir)/manhattan-times-$firstRound.jld", "times")
	newTimes = load("/Users/arthurdelarue/Desktop/VRP-remote/Outputs/$(dir)/manhattan-times-$secondRound.jld", "times")
	man.roadTime = newTimes - oldTimes
	drawManhattan(man, colorScheme = 1)
end

function showAlgorithmResults(dir=DIRECTORY, rounds=1:MAX_ROUNDS)
	man = Manhattan()
	println("**** Original Manhattan ****")
	for i = 1:nv(man.network), j=1:nv(man.network)
		if man.roadTime[i,j] > 0
			man.roadTime[i,j] = 0.038 * man.distances[i,j]
		end
	end
	drawManhattan(man, colorScheme = 0)

	for i in rounds
		println("**** Iteration $i ****")
		man.roadTime = load("/Users/arthurdelarue/Desktop/VRP-remote/Outputs/$(dir)/manhattan-times-$i.jld", "times")
		avg_speed = 0
		num_streets = 0
		for i = 1:nv(man.network), j = out_neighbors(man.network, i)
			avg_speed = (avg_speed * num_streets + man.distances[i,j]/man.roadTime[i,j])/(num_streets + 1)
			num_streets += 1
		end
		avg_speed = 3.6 * avg_speed
		println("Average speed: $avg_speed")
		drawManhattan(man, colorScheme = 0)
	end
end

println("")