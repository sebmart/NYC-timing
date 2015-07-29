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

function print_iis_gurobi(m::Model, im::Gurobi.GurobiMathProgModel)
    """
    Taken from JuMP examples online. Given a LP instance, uses Gurobi to troubleshoot infeasibility.
    Outputs IIS to find "bad" constraints.
    """
    grb = MathProgBase.getrawsolver(im)
    Gurobi.computeIIS(grb)
    numconstr = Gurobi.num_constrs(grb)
    numvar = Gurobi.num_vars(grb)

    iisconstr = Gurobi.get_intattrarray(grb, "IISConstr", 1, numconstr)
    iislb = Gurobi.get_intattrarray(grb, "IISLB", 1, numvar)
    iisub = Gurobi.get_intattrarray(grb, "IISUB", 1, numvar)

    println("Irreducible Inconsistent Subsystem (IIS)")
    println("Variable bounds:")
    lowerBounds = MathProgBase.getvarLB(im)
    upperBounds = MathProgBase.getvarUB(im)
    for i in 1:numvar
        v = Variable(m, i)
        if iislb[i] != 0 && iisub[i] != 0
            println(lowerBounds[i], " <= ", getName(v), " <= ", upperBounds[i])
        elseif iislb[i] != 0
            println(getName(v), " >= ", lowerBounds[i])
        elseif iisub[i] != 0
            println(getName(v), " <= ", upperBounds[i])
        end
    end

    cLowerBounds = MathProgBase.getconstrLB(im)
    cUpperBounds = MathProgBase.getconstrUB(im)
    println("Constraints:")
    for i in 1:numconstr
        if iisconstr[i] != 0
            println(cLowerBounds[i], " <= ", split(string(m.linconstr[i])," = ")[1], " <= ", cUpperBounds[i])
        end
    end
    println("End of IIS")
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
    turnCost::Float64=0.)
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

function findWorstPathIndex(paths::Array{Array{Int}}, numExpensiveTurns::Array{Int}, turnCost::Float64, times::AbstractArray{Float64,2})
    """
    Given some paths (and extra turn costs) figure out which path is the longest
    """
    pathTime = zeros(length(numExpensiveTurns))
    for i = 1:length(pathTime)
        pathTime[i] = sum([times[paths[i][a],paths[i][a+1]] for a = 1:(length(paths[i])-1)]) + numExpensiveTurns[i] * turnCost
    end
    return indmax(pathTime)
end