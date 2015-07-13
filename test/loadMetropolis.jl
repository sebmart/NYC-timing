# loadMetropolis.jl
# Julia 0.4 code to generate metropolis
# Authored by Arthur J Delarue

MAX_VELOCITIES = [1.,1.,1.]#[0.1,0.3,0.5]
MEAN_VELOCITIES = [0.03,0.15,0.3]
OUTSIDE_NODES = [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),(1,8),(2,8),(3,8),(4,8),(5,8),(6,8),(7,8),(8,8),(8,7),(8,6),(8,5),(8,4),(8,3),(8,2),(8,1),(7,1),(6,1),(5,1),(4,1),(3,1),(2,1),(1,1)]
WIDTH = 8
NSUB = 8
TURN_COST = 0.0

using TaxiSimulation
using Distributions
using HDF5, JLD

#Gives the index of vertex, given coords and city number
function coordToLoc(i::Int, j::Int, c::Int, width::Int, nSub::Int)
	subWidth = floor(Int, width/2)
    if c ==0
        return j + (i-1)*width
    else
        return width^2 + (c-1)*subWidth^2 + j + (i-1)*subWidth
    end
end

function load_metropolis_graph(;from_scratch=false)
	path = "/Users/arthurdelarue/.julia/v0.4/TaxiSimulation"
	if from_scratch
		# Create city
		city = Metropolis(WIDTH, NSUB)

		function drawNetworkToText(pb::TaxiProblem, name::String = "graph")
		 	stdin, proc = open(`neato -Tplain -o $(path)/outputs/$(name).txt`, "w")
		 	TaxiSimulation.to_dot(pb,stdin)
		 	close(stdin)
		end

		# GraphViz magic to find coordinates for nodes
		drawNetwork(city)
		GraphViz = Coordinates[]
		indices = Int64[]
		drawNetworkToText(city, "test1")
		fileExists = false
		while (!fileExists)
			sleep(1)
			fileExists = isfile("$(path)/outputs/test1.txt")
		end
		lines = readlines(open("$(path)/outputs/test1.txt"))
		rm("$(path)/outputs/test1.txt")
		index = 2
		while(split(lines[index])[1] == "node")
			push!(indices, convert(Int64, float(split(lines[index])[2])))
			nodeC = Coordinates(float(split(lines[index])[3]), float(split(lines[index])[4]))
			push!(GraphViz, nodeC)
			index += 1
		end

		GraphVizNodes = copy(GraphViz)
		for i = 1:length(GraphViz)
			GraphVizNodes[indices[i]] = GraphViz[i]
		end

		originalNodes = GraphVizNodes
		distances = spzeros(TaxiSimulation.nv(city.network),TaxiSimulation.nv(city.network))

		for i in TaxiSimulation.vertices(city.network), j in TaxiSimulation.out_neighbors(city.network, i)
			distances[i,j] = norm([originalNodes[i].x - originalNodes[j].x, originalNodes[i].y - originalNodes[j].y])
		end

		adjacencyList = Array{Int}[]
		coordinates = zeros(TaxiSimulation.nv(city.network), 2)
		for i = 1:TaxiSimulation.nv(city.network)
			push!(adjacencyList, TaxiSimulation.out_neighbors(city.network, i))
			coordinates[i,1] = originalNodes[i].x
			coordinates[i,2] = originalNodes[i].y
		end

		graph = city.network
		positions = originalNodes

		save("Inputs/input-graph.jld", "adjList", adjacencyList)
		save("Inputs/input-positions.jld", "coordinates", coordinates)
		save("Inputs/input-distances.jld", "distances", distances)
	else
		# load graph
		adjacencyList = load("Inputs/input-graph.jld", "adjList")
		graph = TaxiSimulation.DiGraph(length(adjacencyList))
		for i = 1:TaxiSimulation.nv(graph), j in adjacencyList[i]
			TaxiSimulation.add_edge!(graph, i, j)
		end
		# load coordinates
		coordinates = load("Inputs/input-positions.jld", "coordinates")
		positions = TaxiSimulation.Coordinates[]
		for i = 1:TaxiSimulation.nv(graph)
			push!(positions, TaxiSimulation.Coordinates(coordinates[i,1], coordinates[i,2]))
		end
		# load distances
		distances = load("Inputs/input-distances.jld", "distances")
	end
	return graph, positions, distances
end

function define_travel_times(;from_scratch = false, turn_cost::Float64 = TURN_COST)
	input_speed = 0
	function gen_speed(index::Int)
		input_speed += 0.001 #MEAN_VELOCITIES[index]
		speed = -1
		while speed < 0 || speed > 9.9#MAX_VELOCITIES[index]
			speed = rand(Normal(input_speed, input_speed/3))#15))
		end
		return input_speed
	end
	if from_scratch
		graph, positions, distances = load_metropolis_graph(from_scratch=false)

		minRoadTime = spzeros(TaxiSimulation.nv(graph),TaxiSimulation.nv(graph))
		meanRoadTime = spzeros(TaxiSimulation.nv(graph),TaxiSimulation.nv(graph))
		# Set normal roads
		for i in TaxiSimulation.vertices(graph), j in TaxiSimulation.out_neighbors(graph, i)
			minRoadTime[i,j] = distances[i,j]/MAX_VELOCITIES[2]
			meanRoadTime[i,j] = distances[i,j]/gen_speed(2)
		end
		# Set fast roads on outside of main city
		for k = 1:(length(OUTSIDE_NODES)-1)
			i1,j1 = OUTSIDE_NODES[k]
			node1 = coordToLoc(i1,j1,0,WIDTH,NSUB)
			i2,j2 = OUTSIDE_NODES[k+1]
			node2 = coordToLoc(i2,j2,0,WIDTH,NSUB)
			minRoadTime[node1,node2] = distances[node1,node2]/MAX_VELOCITIES[3]
			minRoadTime[node2,node1] = distances[node2,node1]/MAX_VELOCITIES[3]
			meanRoadTime[node1,node2] = distances[node1,node2]/gen_speed(3)
			meanRoadTime[node2,node1] = distances[node2,node1]/gen_speed(3)
		end
		# Set slow roads between suburbs and city
		for k = 1:8
			suburbNode = coordToLoc(1,1,k,WIDTH,NSUB)
			cityNode = minimum(TaxiSimulation.out_neighbors(graph, suburbNode))
			minRoadTime[suburbNode, cityNode] = distances[suburbNode, cityNode]/MAX_VELOCITIES[1]
			minRoadTime[cityNode, suburbNode] = distances[cityNode, suburbNode]/MAX_VELOCITIES[1]
			meanRoadTime[suburbNode, cityNode] = distances[suburbNode, cityNode]/gen_speed(1)
			meanRoadTime[cityNode, suburbNode] = distances[cityNode, suburbNode]/gen_speed(1)
		end

		cityPaths = realPaths(graph, meanRoadTime, meanRoadTime, positions, turn_cost, turn_cost)
		traveltime = cityPaths.traveltime

		save("Inputs/input-realTimes.jld", "realTimes", cityPaths.traveltime)
		save("Inputs/input-speedLimits.jld", "speedLimits", minRoadTime)
		save("Inputs/input-meanTimes.jld", "meanTimes", meanRoadTime)
	else
		traveltime = load("Inputs/input-realTimes.jld", "realTimes")
		minRoadTime = load("Inputs/input-speedLimits.jld", "speedLimits")
		meanRoadTime = load("Inputs/input-meanTimes.jld", "meanTimes")
	end
	return traveltime, minRoadTime, meanRoadTime
end

function generate_rides(prob = 0.7)
	graph, positions, distances = load_metropolis_graph(from_scratch=false)
	traveltime, minRoadTime, meanRoadTime = define_travel_times(from_scratch=false)

	out = [copy(TaxiSimulation.out_neighbors(graph, i)) for i = TaxiSimulation.vertices(graph)]
	out2 = [vcat([copy(TaxiSimulation.out_neighbors(graph, j)) for j in out[i]]) for i = TaxiSimulation.vertices(graph)]

	cityToCityDist = Categorical([prob, (1-prob)/6, (1-prob)/6, (1-prob)/6, (1-prob)/6, (1-prob)/6, (1-prob)/6])
	cityToSuburbDist = Categorical([prob, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10])
	suburbToCityDist = Categorical([prob, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10, (1-prob)/10])
	suburbToSuburbDist = Categorical([prob, (1-prob)/4, (1-prob)/4, (1-prob)/4, (1-prob)/4])

	numRides = zeros(Int, (TaxiSimulation.nv(graph), TaxiSimulation.nv(graph)))
	travelTimes = zeros(TaxiSimulation.nv(graph), TaxiSimulation.nv(graph))
	for i = 1:TaxiSimulation.nv(graph), j = 1:TaxiSimulation.nv(graph)
		if i != j && !(j in out[i]) && !(j in out2[i])
			# city to city
			if i < 65 && j < 65
				numRides[i,j] = rand(cityToCityDist) - 1
			# city to suburb
			elseif i < 65 && j >= 65
				numRides[i,j] = rand(cityToSuburbDist) - 1
			# suburb to city
			elseif i >= 65 && j < 65
				numRides[i,j] = rand(suburbToCityDist) - 1
			# suburb to suburb
			else
				numRides[i,j] = rand(suburbToSuburbDist) - 1
			end
			if numRides[i,j] > 0
				travelTimes[i,j] = mean(rand(Normal(traveltime[i,j], sqrt(traveltime[i,j])), numRides[i,j]))
				while travelTimes[i,j] < 0
					travelTimes[i,j] = mean(rand(Normal(traveltime[i,j], sqrt(traveltime[i,j])), numRides[i,j]))
				end
			end
		end
	end

	save("Inputs/input-travelTimes-$(prob).jld", "travelTimes", travelTimes)
	save("Inputs/input-numRides-$(prob).jld", "numRides", numRides)
	return travelTimes, numRides
end

function create_metropolis(graph::TaxiSimulation.DiGraph, roadTime::AbstractArray{Float64,2})
	city = Metropolis(WIDTH, NSUB)
	city.network = graph
	city.roadTime = roadTime
	return city
end

function updatedMetropolis(prob::Float64; max_rounds = 20, num_iter = -1)
	graph, positions, distances = load_metropolis_graph()
	if num_iter < 0 
		if isfile("Outputs/minr1_i$(max_rounds)_metropolis_$(prob)/end.jld")
			num_iter = load("Outputs/minr1_i$(max_rounds)_metropolis_$(prob)/end.jld", "num_iter")
		else
			num_iter = max_rounds
		end
	end
	roadTime = load("Outputs/minr1_i$(max_rounds)_metropolis_$(prob)/metropolis-times-$(num_iter).jld", "times")
	return create_metropolis(graph, roadTime)
end

function showRealMetropolis()
	g, p, d = load_metropolis_graph()
	tt, mt, at = define_travel_times()
	visualize(create_metropolis(g, at))
end

function showStartingMetropolis()
	g, p, d = load_metropolis_graph()
	tt, mt, at = define_travel_times()
	visualize(create_metropolis(g, mt))
end