#----------------------------------------
#-- Set all the types
#----------------------------------------


"""
to keep?
"""
immutable ShortPaths
  traveltime::Array{Float64,2}
  travelcost::Array{Float64,2}
  previous::Array{Int,2}
end

ShortPaths() = ShortPaths( Array(Float64, (0,0)), Array(Float64, (0,0)), Array(Int, (0,0)))

immutable RealPaths
  previous::Array{Int,2}
  traveltime::Array{Float64,2}
  real_destinations::Array{Int, 2}
end

# define heap entry data type
immutable DijkstraEntry{Float64}
  vertex::Int
  dist::Float64
  cost::Float64
end

# define appropriate comparators for heap entries
<(e1::DijkstraEntry, e2::DijkstraEntry) = e1.dist < e2.dist
Base.isless(e1::DijkstraEntry, e2::DijkstraEntry) = e1.dist < e2.dist
