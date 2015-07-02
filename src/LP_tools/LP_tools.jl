# LP_tools.jl
# Contains helpful functions to debug Gurobi output
# Authored By Arthur J Delarue on 6/10/15

function print_iis_gurobi(m::Model)
	"""
	Taken from JuMP examples online. Given a LP instance, uses Gurobi to troubleshoot infeasibility.
	Outputs IIS to find "bad" constraints.
	"""
    grb = MathProgBase.getrawsolver(getInternalModel(m))
    Gurobi.computeIIS(grb)
    numconstr = Gurobi.num_constrs(grb)
    numvar = Gurobi.num_vars(grb)

    iisconstr = Gurobi.get_intattrarray(grb, "IISConstr", 1, numconstr)
    iislb = Gurobi.get_intattrarray(grb, "IISLB", 1, numvar)
    iisub = Gurobi.get_intattrarray(grb, "IISUB", 1, numvar)

    println("Irreducible Inconsistent Subsystem (IIS)")
    println("Variable bounds:")
    for i in 1:numvar
        v = Variable(m, i)
        if iislb[i] != 0 && iisub[i] != 0
            println(getLower(v), " <= ", getName(v), " <= ", getUpper(v))
        elseif iislb[i] != 0
            println(getName(v), " >= ", getLower(v))
        elseif iisub[i] != 0
            println(getName(v), " <= ", getUpper(v))
        end
    end

    println("Constraints:")
    for i in 1:numconstr
        if iisconstr[i] != 0
            println(m.linconstr[i])
        end
    end
    println("End of IIS")
end

function parallelShortestPaths(n::Network, roadTime::SparseMatrixCSC{Float64, Int},roadCost::SparseMatrixCSC{Float64, Int})
	"""
	Computes all-pairs shortest paths in a network by distributing jobs evenly over avilable cores. Not fastest implementation because jobs are divided into as many pieces as there are CPUs, and each CPU ends up with a large chunk of jobs, so the time is limited by the slowest processor.
	"""
	# Get number of nodes
	nLocs  = length( vertices(n))
	# Find number of cores
	nCpus = length(workers())

	# Start processes
	results = cell(nCpus)
	for i = 1:nCpus
		results[i] = @spawnat i partialShortestPaths(n, (div((i-1) * nLocs, nCpus)+1):(div(i * nLocs,nCpus)), roadTime, roadCost)
	end
	# Create return arrays
	pathTime = Array(Float64, (nLocs,nLocs))
	pathCost = Array(Float64, (nLocs,nLocs))
	previous = Array(Int, (nLocs,nLocs))
	# Fetch results and put them together in big array
	for i = 1:length(results)
		tmp_time, tmp_cost, tmp_prev = fetch(results[i])
		irange = (div((i-1) * nLocs, nCpus)+1):(div(i * nLocs,nCpus))
		pathTime[irange,:] = tmp_time
		pathCost[irange,:] = tmp_cost
		previous[irange,:] = tmp_prev
	end
	
	return ShortPaths(pathTime, pathCost, previous)
end

@everywhere function parallelShortestPathsAuto(n::Network, roadTime::SparseMatrixCSC{Float64, Int},roadCost::SparseMatrixCSC{Float64, Int}; batch_size::Int = 400)
	"""
	Computes all-pairs shortest paths in a network by separating jobs into chunks of size batch_size (default 100) and letting pmap() deal with assigning them to the next available CPU. With the right batch_size this is the fastest method.
	"""
	# Get number of nodes
	nLocs  = length( vertices(n))
	# Initialize arrays for shortest paths
	pathTime = Array(Float64, (nLocs,nLocs))
	pathCost = Array(Float64, (nLocs,nLocs))
	previous = Array(Int, (nLocs,nLocs))

	# Split locations into chunks that can be parallelized
	locations = splitLocations(1:nLocs, batch_size)

	# Wrapper function, created so that we only need to pass one parameter to pmap
	function multiDijkstraInstance(nodes::UnitRange{Int})
		return partialShortestPaths(n, nodes, roadTime, roadCost)
	end
	# Call pmap on wrapper function
	results = pmap(multiDijkstraInstance, locations)

	# Process output of pmap and return
	for i in 1:length(results)
		tmp_time, tmp_cost, tmp_prev = results[i]
		irange = locations[i]
		pathTime[irange,:] = tmp_time
		pathCost[irange,:] = tmp_cost
		previous[irange,:] = tmp_prev
	end
	return ShortPaths(pathTime, pathCost, previous)
end

@everywhere function parallelShortestPathsWithTurnsAuto(n::Network, new_n::Network, newRoadTime::SparseMatrixCSC{Float64, Int},new_nodes::Array{Array{Int}}; batch_size::Int = 400)
    """
    Computes all-pairs shortest paths with turns in a properly modified network by separating jobs into chunks of size batch_size (default 100) and letting pmap() deal with assigning them to the next available CPU. With the right batch_size this is the fastest method.
    """
    # Get number of nodes
    nLocs  = nv(n)
    nRealLocs = nv(new_n)
    # Initialize arrays for shortest paths
    pathTime = Array(Float64, (nLocs,nLocs))
    realDsts = Array(Int, (nLocs,nLocs))
    previous = Array(Int, (nLocs,nRealLocs))

    # Split locations into chunks that can be parallelized
    locations = splitLocations(1:nLocs, batch_size)

    # Wrapper function, created so that we only need to pass one parameter to pmap
    function multiDijkstraInstance(nodes::UnitRange{Int})
        return partialShortestPathsWithTurns(n, new_n, nodes, newRoadTime, new_nodes)
    end
    # Call pmap on wrapper function
    results = pmap(multiDijkstraInstance, locations)

    # Process output of pmap and return
    for i in 1:length(results)
        tmp_dsts, tmp_time, tmp_prev = results[i]
        irange = locations[i]
        pathTime[irange,:] = tmp_time
        realDsts[irange,:] = tmp_dsts
        previous[irange,:] = tmp_prev
    end
    return RealPaths(previous, pathTime, realDsts)
end

@everywhere function splitLocations(locs::UnitRange{Int}, size::Int)
	"""
	Function divides a unit range into an array of smaller unitranges of a given size, such that all the same integers are still covered.
	"""
    result = UnitRange{Int}[]
    b = 1
    # While there are still elements to be covered
    while b <= length(locs)
    	# Deal with the fact that the last unit range may slightly shorter
        finish = min(b+size-1, length(locs))
      	# Append to list
        push!(result, b:finish)
        b += size
    end
    return result
end

@everywhere function partialShortestPaths(n::Network, locs::UnitRange{Int}, roadTime::SparseMatrixCSC{Float64, Int}, roadCost::SparseMatrixCSC{Float64, Int})
	"""
	Given a graph and a set of locations in the graph, runs Dijkstra from each provided node to find the SSSP from each given node (subset of APSP).
	"""
	# Get number of nodes
	nLocs  = length( vertices(n))
	# Initialize return types
	partialPathTime = zeros(Float64, (length(locs),nLocs))
	partialPathCost = zeros(Float64, (length(locs),nLocs))
	partialPrevious = zeros(Int, (length(locs),nLocs))

	# Run dijkstra from each node
	for i in locs
		parents, dists, costs = custom_dijkstra_par(n, i, roadTime, roadCost)
		partialPrevious[i-minimum(locs)+1,:] = parents
		partialPathTime[i-minimum(locs)+1,:] = dists
		partialPathCost[i-minimum(locs)+1,:] = costs
	end
	# Output result
	return partialPathTime, partialPathCost, partialPrevious
end

@everywhere function partialShortestPathsWithTurns(n::Network, new_n::Network, locs::UnitRange{Int}, newRoadTime::SparseMatrixCSC{Float64,Int}, new_nodes::Array{Array{Int}})
    """
    Given a graph, a modified graph (for turn costs) and a set of locations in the graph, runs Dijkstra from each provided node to find the SSSP from each given node (subset of APSP).
    """
    # Get number of nodes
    nLocs  = nv(n)
    nRealLocs = nv(new_n)
    # Initialize return types
    partialPathTime = zeros(Float64, (length(locs),nLocs))
    partialRealDsts = zeros(Int, (length(locs),nLocs))
    partialPrevious = zeros(Int, (length(locs),nRealLocs))

    # Run dijkstra from each node
    for i in locs
        parents, dists, real_dsts = custom_dijkstra_with_turn_cost(n, new_n, i, newRoadTime, new_nodes)
        partialPrevious[i-minimum(locs)+1,:] = parents
        partialPathTime[i-minimum(locs)+1,:] = dists
        partialRealDsts[i-minimum(locs)+1,:] = real_dsts
    end
    # Output result
    return partialRealDsts, partialPathTime, partialPrevious
end

@everywhere function custom_dijkstra_par(
    g::AbstractGraph,
    src::Int,
    edge_dists::AbstractArray{Float64, 2},
    edge_costs::AbstractArray{Float64, 2})
    """
    Performs Dijkstra's algorithm from node src on graph g using weights edge_dists.
    Returns parents, travel times and travel costs as vectors.
    """
    nvg = nv(g)
    dists = fill(typemax(Float64), nvg)
    costs = fill(typemax(Float64), nvg)
    parents = zeros(Int, nvg)
    visited = falses(nvg)
    h = DijkstraEntry{Float64}[]
    sizehint(h, nvg)
    H = mutable_binary_minheap(h)
    hmap = zeros(Int, nvg)
    dists[src] = 0.0
    costs[src] = 0.0
    ref = push!(H, DijkstraEntry{Float64}(src, dists[src], costs[src]))
    hmap[src] = ref
    while !isempty(H)
        hentry = pop!(H)
        u = hentry.vertex
        for v in out_neighbors(g,u)
            if dists[u] == typemax(Float64)
                alt = typemax(Float64)
                alt2 = typemax(Float64)
            else
                alt = dists[u] + edge_dists[u,v]
                alt2 = costs[u] + edge_costs[u,v]
            end
            if !visited[v]
                dists[v] = alt
                costs[v] = alt2
                parents[v] = u
                visited[v] = true
                ref = push!(H, DijkstraEntry{Float64}(v, alt, alt2))
                hmap[v] = ref
            else
                if alt < dists[v]
                    dists[v] = alt
                    costs[v] = alt2
                    parents[v] = u
                    update!(H, hmap[v], DijkstraEntry{Float64}(v, alt, alt2))
                end
            end
        end
    end

    dists[src] = 0
    costs[src] = 0.0
    parents[src] = 0

    return parents, dists, costs
end

@everywhere function custom_dijkstra_with_turn_cost(
    g::AbstractGraph,
    newGraph::AbstractGraph,
    src::Int,
    new_edge_dists::AbstractArray{Float64, 2},
    new_nodes::Array{Array{Int}})
    """
    Modified version of Dijkstra, to be run on a modified graph (that includes left turn costs). Takes in old graph g, as well as modified graph newGraph. Runs Dijkstra from source src, using edge weights new_edge_dists.
    new_nodes : map from nodes in g to nodes in newGraph
    Returns:
    parents, a vector of parents in the MODIFIED graph
    real_dists/real_costs, vectors of distances/costs in the UNMODIFIED graph
    real_dsts, a map of each node to the corresponding reached node in Dijkstra, in the MODIFIED graph
    """
    n = nv(g)
    real_dists = zeros(Float64, n)
    real_dsts = zeros(Int, n)

    # Compute shortest paths in new graph
    nvg = nv(newGraph)
    dists = fill(typemax(Float64), nvg)
    parents = zeros(Int, nvg)
    visited = falses(nvg)
    # Create heap
    h = DijkstraEntry{Float64}[]
    sizehint(h, nvg)
    H = mutable_binary_minheap(h)
    hmap = zeros(Int, nvg)
    
    # Add all sources to heap (possibly multiple)
    for src2 in new_nodes[src]
        dists[src2] = 0.0
        ref = push!(H, DijkstraEntry{Float64}(src2, dists[src2], 0.0))
        hmap[src2] = ref
        visited[src2] = true
    end

    # Run Dijkstra normally
    while !isempty(H)
        hentry = pop!(H)
        u = hentry.vertex
        for v in out_neighbors(newGraph,u)
            if dists[u] == typemax(Float64)
                alt = typemax(Float64)
            else
                alt = dists[u] + new_edge_dists[u,v]
            end
            if !visited[v]
                dists[v] = alt
                parents[v] = u
                visited[v] = true
                ref = push!(H, DijkstraEntry{Float64}(v, alt, 0.0))
                hmap[v] = ref
            else
                if alt < dists[v]
                    dists[v] = alt
                    parents[v] = u
                    update!(H, hmap[v], DijkstraEntry{Float64}(v, alt, 0.0))
                end
            end
        end
    end

    # Update back to original graph
    for node = 1:n
        real_dists[node], index = findmin(dists[new_nodes[node]])
        real_dsts[node] = index + new_nodes[node][1] - 1
    end

    return parents, real_dists, real_dsts
end

@everywhere function modifyGraphForDijkstra(
    g::AbstractGraph,
    edge_dists::AbstractArray{Float64, 2},
    coords::Array{Coordinates, 1};
    turn_cost::Float64 = 10.0)
    """
    Given a graph g, modifies it so that left turns are afflicted with the extra cost turn_cost.
    Requires old edge_costs as well as geographic coordinates of nodes (to determine left turns).
    Returns new network, new edge weights, and map from nodes in g to nodes in the new graph as a list of lists.
    """

    function getAngleToPoint(currentAngle::Float64, currentX::Float64, currentY::Float64, targetX::Float64, targetY::Float64)
        """
        Get angle to a given point from a given point, given a current heading.
        """
        angle = atan2(targetY - currentY, targetX - currentX) - currentAngle
        # Make sure it's between -pi and pi
        while angle > pi
            angle = angle - 2 * pi
        end
        while angle < -pi
            angle = angle + 2 * pi
        end
        return angle
    end

    function getAngleEdges(xs::Float64, ys::Float64, xm::Float64, ym::Float64, xe::Float64, ye::Float64)
        """
        Get angle between two edges in a graph.
        """
        currentAngle = getAngleToPoint(0.0, xs, ys, xm, ym)
        edgeAngle = getAngleToPoint(currentAngle, xm, ym, xe, ye)
        return edgeAngle
    end
    
    # Define some early variables
    nvg = nv(g)
    out = [sort(out_neighbors(g,i)) for i = 1:nvg]
    inn = [sort(in_neighbors(g,i)) for i = 1:nvg]
    # Create new graph
    newGraph = Network()
    new_nodes = Array{Int}[]
    sizehint(new_nodes, nvg)
    # Deal with vertices first
    for i in 1:nvg
        push!(new_nodes,Int[])
        sizehint(new_nodes[i], length(inn[i]))
        # Create as many new vertices as there are edges into the original vertex
        for j = 1:length(inn[i])
            id = add_vertex!(newGraph)
            push!(new_nodes[i], id)
        end
    end
    new_edge_dists = spzeros(nv(newGraph), nv(newGraph))
    # Now deal with edges
    for i = 1:nvg, j = 1:length(inn[i]), k = 1:length(out[i])
        # Find correct sub-node of i to connect to k
        l = findfirst(inn[out[i][k]], i)
        add_edge!(newGraph, new_nodes[i][j], new_nodes[out[i][k]][l])
        src = inn[i][j]
        dst = out[i][k]
        # Add extra edge weight appropriately
        # Make sure not to add weight if the node has only one incoming edge and one outgoing edge
        angle = getAngleEdges(coords[src].x, coords[src].y, coords[i].x, coords[i].y, coords[dst].x, coords[dst].y)
        if angle < pi/4 || (length(out[i]) == 1 && length(inn[i]) == 1)
            new_edge_dists[new_nodes[i][j], new_nodes[dst][l]] = edge_dists[i, dst]
        else
            new_edge_dists[new_nodes[i][j], new_nodes[dst][l]] = edge_dists[i, dst] + turn_cost
        end
    end
    return newGraph, new_edge_dists, new_nodes
end

@everywhere function getInverseMapping(new_nodes::Array{Array{Int}}, numNewVertices::Int)
    """
    Given mapping from old nodes to new, returns mapping from new nodes to old.
    """
    old_nodes = zeros(Int, numNewVertices)
    for (i, array) in enumerate(new_nodes)
        for j in array
            old_nodes[j] = i
        end
    end
    return old_nodes
end

@everywhere function splitArray(array::Array{Int}, batch_size::Int)
    """
    Given an array, returns an array of arrays of size batch_size containing all the same elements.
    """
    indices = splitLocations(1:length(array), batch_size)
    split_array = fill(Int[],length(indices))
    for i = 1:length(indices)
        split_array[i] = copy(array[indices[i]])
    end
    return split_array
end

@everywhere function reconstructMultiplePathsWithExpensiveTurnsParallel(
    previous::Array{Int, 2},
    srcs::Array{Int},
    dsts::Array{Int},
    old_nodes::Array{Int},
    real_dsts::Array{Int},
    old_edge_dists::AbstractArray{Float64,2},
    new_edge_dists::AbstractArray{Float64,2};
    batch_size=30000)
    """
    Given an array of sources and a corresponding array of destinations, returns the used paths in the old graph, and the number of expensive turns (i.e. left turns) along each path.
    Parallelized for speed.
    """
    paths = fill(Int[], length(srcs))
    expensiveTurns = zeros(Int, length(srcs))

    # Split array of sources
    segmentedSources = splitArray(srcs, batch_size)
    segmentedDests = splitArray(dsts, batch_size)
    indices = splitLocations(1:length(srcs), batch_size)

    function callReconstruction(orig::Array{Int}, dest::Array{Int})
        return reconstructMultiplePathsWithExpensiveTurns(previous, orig, dest, old_nodes, real_dsts, old_edge_dists, new_edge_dists)
    end

    # Call pmap on wrapper function
    results = pmap(callReconstruction, segmentedSources, segmentedDests)

    # Process output of pmap and return
    for i = 1:length(results)
        tmp_paths, tmp_turns = results[i]
        paths[indices[i]] = tmp_paths
        expensiveTurns[indices[i]] = tmp_turns
    end

    return paths, expensiveTurns
end

@everywhere function reconstructMultiplePathsWithExpensiveTurns(
    previous::Array{Int,2},
    srcs::Array{Int},
    dsts::Array{Int},
    old_nodes::Array{Int},
    real_dsts::Array{Int,2},
    old_edge_dists::AbstractArray{Float64,2},
    new_edge_dists::AbstractArray{Float64,2})
    """
    Given an array of sources and a corresponding array of destinations, returns the used paths in the old graph, and the number of expensive turns (i.e. left turns) along each path.
    """
    expensiveTurns = zeros(Int, length(srcs))
    paths = fill(Int[], length(srcs))
    for i = 1:length(srcs)
        tmp_path, tmp_expensiveTurns = reconstructPathWithExpensiveTurns(previous[srcs[i],:], srcs[i], dsts[i], old_nodes, real_dsts[srcs[i],:], old_edge_dists, new_edge_dists)
        paths[i] = tmp_path
        expensiveTurns[i] = tmp_expensiveTurns
    end
    return paths, expensiveTurns
end

@everywhere function reconstructPathWithExpensiveTurns(
    previous::Array{Int},
    src::Int,
    dst::Int,
    old_nodes::Array{Int},
    real_dsts::Array{Int},
    old_edge_dists::AbstractArray{Float64,2},
    new_edge_dists::AbstractArray{Float64,2})
    """
    Given a previous object, a source and a destination, returns the used path in the old graph, and the number of expensive turns (i.e. left turns) along the path.
    """
    expensiveTurns = 0::Int
    # Initialize path in real graph
    pathNodes = [dst]
    k = real_dsts[dst]
    # Initialize path in modified graph (only keep last vertex)
    previousActualPathNode = k
    sizehint(pathNodes, 200)
    # Look for path in previous matrix
    while old_nodes[k] != src
        k = previous[k]
        # Check if edge was modified (i.e. by presence of left turn)
        if new_edge_dists[k, previousActualPathNode] - old_edge_dists[old_nodes[k], pathNodes[end]] > 1e-7
            expensiveTurns += 1
        end
        # Add node to path
        push!(pathNodes, old_nodes[k])
        # Update previous actual node
        previousActualPathNode = k
    end
    reverse!(pathNodes)
    return pathNodes, expensiveTurns
end

function reconstructFullPath(
    previous::Array{Int},
    src::Int,
    dst::Int,
    old_nodes::Array{Int},
    real_dsts::Array{Int},
    old_edge_dists::AbstractArray{Float64,2},
    new_edge_dists::AbstractArray{Float64,2})
    """
    Return full path as list of (road, time) tuples. 
    """
    # Initialize path in real graph
    path = (Road, Float64)[]
    k = real_dsts[dst]
    # Initialize path in modified graph (only keep last vertex)
    previousActualPathNode = k
    previousPathNode = dst
    # Look for path in previous matrix
    while old_nodes[k] != src
        k = previous[k]
        # Add edge to path
        push!(path, (Road(old_nodes[k],previousPathNode),old_edge_dists[old_nodes[k], previousPathNode]))
        # Check if edge was modified (i.e. by presence of left turn)
        if new_edge_dists[k, previousActualPathNode] - old_edge_dists[old_nodes[k], previousPathNode] > 1e-7
            push!(path, (Road(old_nodes[k],old_nodes[k]), new_edge_dists[k, previousActualPathNode] - old_edge_dists[old_nodes[k], previousPathNode]))
        end
        # Update previous actual node
        previousActualPathNode = k
        previousPathNode = old_nodes[k]
    end
    reverse!(path)
    return path
end

function writeDataToFile(directory::String,
    model_type::String="",
    max_rounds::Int=0,
    min_rides::Int=0,
    radius::Int=0,
    regularization::Bool=false,
    preprocess::Bool=false,
    num_clusters::Int=0,
    sample_size::Int=0,
    turnCost::Float64=0.,
    lambda::Float64=0.)
    """
    Writes relevant information for the LP to a file in the provided directory.
    """
    fileName = string(directory, "/info.txt")
    outputFile = open(fileName, "w")
    write(outputFile, string(fileName, "\n"))
    write(outputFile, string(now(), "\n"))
    write(outputFile, string("Model type    :\t", model_type, "\n"))
    write(outputFile, string("Max iterations:\t", max_rounds, "\n"))
    write(outputFile, string("Min # of rides:\t", min_rides, "\n"))
    write(outputFile, string("Node radius   :\t", radius, "m\n"))
    if regularization
        write(outputFile, "Regularization:\tON\n")
        write(outputFile, string("Lambda        :\t",lambda, "\n"))
    else
        write(outputFile, "Regularization:\tOFF\n")
        write(outputFile, string("Lambda        :\tN/A\n"))
    end
    if preprocess
        write(outputFile, "Preprocessing :\tON\n")
        write(outputFile, string("# of clusters :\t", num_clusters, "\n"))
        write(outputFile, string("Max samp. size:\t", sample_size, "\n"))
    else
        write(outputFile, "Preprocessing :\tOFF\n")
        write(outputFile, "# of clusters :\tN/A\n")
        write(outputFile, "Max samp. size:\tN/A\n")
    end
    write(outputFile, string("Left turn cost:\t", turnCost, "s\n"))
    close(outputFile)
end

function writeDataToFile(directory::String,
    model_type::String="",
    max_rounds::Int=0,
    min_rides::Int=0,
    radius::Int=0,
    preprocess::Bool=false,
    num_clusters::Int=0,
    sample_size::Int=0,
    turnCost::Float64=0.,
    lambda::Float64=0.)
    """
    Writes relevant information for the LP to a file in the provided directory.
    """
    fileName = string(directory, "/info.txt")
    outputFile = open(fileName, "w")
    write(outputFile, string(fileName, "\n"))
    write(outputFile, string(now(), "\n"))
    write(outputFile, string("Model type    :\t", model_type, "\n"))
    write(outputFile, string("Max iterations:\t", max_rounds, "\n"))
    write(outputFile, string("Min # of rides:\t", min_rides, "\n"))
    write(outputFile, string("Node radius   :\t", radius, "m\n"))
    write(outputFile, string("Lambda        :\t",lambda, "\n"))
    if preprocess
        write(outputFile, "Preprocessing :\tON\n")
        write(outputFile, string("# of clusters :\t", num_clusters, "\n"))
        write(outputFile, string("Max samp. size:\t", sample_size, "\n"))
    else
        write(outputFile, "Preprocessing :\tOFF\n")
        write(outputFile, "# of clusters :\tN/A\n")
        write(outputFile, "Max samp. size:\tN/A\n")
    end
    write(outputFile, string("Left turn cost:\t", turnCost, "s\n"))
    close(outputFile)
end

function appendDataToFile(directory::String,
    numBrokenConstraints::Int)
    """
    Adds the number of broken constraints through the relaxation parameter to info.txt
    """
    fileName = string(directory, "/info.txt")
    outputFile = open(fileName, "a")
    write(outputFile, string("Constraints violated by regularization: ", numBrokenConstraints, "\n"))
    close(outputFile)
end
