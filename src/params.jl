# params.jl
# Select parameters and save them into JSON files
# Authored by Arthur Delarue on 10/28/15

using JSON

PARAMS = Dict(
"MODEL" => "relaxed_smooth",			# relaxed/strict _ smooth/nothing/random
"LAST_SMOOTH" => false, 				# true if last iteration is smooth, false otherwise. Ignore if model is smooth already
"ERROR_COMPUTATION" => "both", 			# average/single-ride/both

"MAX_ROUNDS" => 50,						# max number of iterations

"MIN_RIDES" => 1,						# min number of rides
"RADIUS" => 200,						# radius (4D Euclidian)
"YEAR" => 2013,						
"START_TIME" => 12,
"END_TIME" => 14,
"START_MONTH" => 1,
"END_MONTH" => 12,
"WEEKDAYS" => true,						# true for weekdays, false for weekends
"METHOD" => "average",					# average/nearest [neighbor]

"RANDOM_CONSTRAINTS" => false, 			# true for purely random constraints, false otw
"DYNAMIC_CONSTRAINTS" => true, 			# true for dynamic constraints, false otw
"GLOBAL_CONSTRAINT_UPDATE" => true,		# true if we remove worst paths globally, false otw

"SAMPLE_SIZE" => 1000, 					# starting number of constraints
"NUM_OD_ADDED" => 2500,					# number of (O,D) pairs to add
"NUM_OD_REMOVED" => 1000,				# number of (O,D) pairs to remove
"UPDATE_EVERY_N_ITERATIONS" => 2, 		# add number of (O,D) above every $N iterations
"ADD_IF_ABOVE_PERCENTILE" => 0.7,		# only add rides that have an error above this percentile
"REMOVE_IF_BELOW_PERCENTILE" => 0.3,	# only remove rides that have an error below this percentile

"TURN_COST" => 0.,	 					# turning cost initial value
"TURN_COST_AS_VARIABLE" => false, 		# true if LP updates turning cost, false otw

"DELTA_BOUND" => [0.06, 0.05,0.04],		# relaxation base (sketchy)

"MAX_NUM_PATHS_PER_OD" => 3, 			# at least 1 please

"COMPUTE_FINAL_SHORTEST_PATHS" => true, # Should always be true unless trying to save time somehow

"START_SIMPLE" => false, 				# true if initial simple LP is used, false otw

"METROPOLIS" => false 					# always set this to false unless you wish to crash everything
)

function create_JSON_file(params::Dict{ASCIIString,Any}, fileName::AbstractString)
	file = open(fileName, "w")
	write(file, json(params))
	close(file)
end

function add_JSON_file(params::Dict{ASCIIString,Any})
	i = 0
	while isfile("JSONparams/params-$i.json")
		i += 1
	end
	create_JSON_file(params, "JSONparams/params-$i.json")
end

add_JSON_file(PARAMS)