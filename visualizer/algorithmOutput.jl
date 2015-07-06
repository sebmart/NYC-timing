# algorithmOutput.jl
# Use manhattanVisualizations to see step-by-step progression of algorithm
# Authored by Arthur J Delarue on 6/16/15

using SFML, LightGraphs, HDF5, JLD
using Base.Collections

include("../src/LP_tools/definitions.jl")
include("../src/LP_tools/manhattan.jl")

include("manhattanVisualization.jl")

DIRECTORY = "new_r40_minr3_i3_wd_1214_quad_clust50_rides5000"
MAX_ROUNDS = 3

function showAlgorithmResults(dir=DIRECTORY, rounds=1:MAX_ROUNDS)
	man = load("../src/Cities/Saved/full_manhattan.jld","pb")
	println("**** Original Manhattan ****")
	drawManhattan(man)

	for i = rounds
		println("**** Iteration $i ****")
		man.roadTime = load("../src/Outputs/$(dir)/manhattan-times-$i.jld", "times")
		drawManhattan(man)
	end
end

println("")