# NYC-timing

Code for travel time estimation from taxi ride data.

## Workflow
#### Interactive mode
- Start julia with required number of cores (usually 8 or more): 
```bash
$ julia -p 8
```
- Load all dependencies
```julia
julia> include("intro.jl")
```
- Run LP method
```julia
julia> include("fast_LP.jl")
```
#### Using qsub
- Load parameter files that you want to run by first editing src/params.jl, then running
```bash
julia params.jl
```
This creates a parameter file in the JSONparams/ directory.
- Choose desired parameter files by listing them in JSONparams/paramList.txt
- Submit jobs
```bash
./script.sh
```

## File index
Files not listed here but present in repository are not useful.
#### Main files
- src/intro.jl: Loads all dependencies and helper functions into memory
- src/fast_LP.jl: Run LP methods (more detailed description below)
- src/simple_LP.jl: Simple one-variable LP (to be used as initial step for main method)
- src/call_fast_LP.jl: Self-contained LP calling file, called by bash script to be fed directly into computing grid. Prameters are read from JSON file.
- src/juliaCall.sh: bash subroutine that calls call_fast_LP.jl
- src/script.sh: reads list of parameter files, submits a job for each one
- src/params.jl: creates JSON file for desired parameters

#### Helper files
- src/LP_tools/LP_tools.jl: miscellaneous helper functions (especially saving output and error info to files)
- src/LP_tools/data_handling.jl: most recent loading of constraints using Sébastien's dataframe method
- src/LP_tools/definitions.jl: copied from TaxiSimulation (with a couple extra things)
- src/LP_tools/dynamic_constraints.jl: helper functions for dynamic (O,D) addition (and deletion... pending)
- src/LP_tools/load_constraints.jl: contained functions to load Manhattan and taxi data into memory. Now only used to load Manhattan (all data loading handled by data_handling.jl)
- src/LP_tools/manhattan.jl: copied from TaxiSimulation with a couple of additions (like ability to create an empty Manhattan object)
- src/LP_tools/parallel_tools.jl: helper functions for parallelized operations (shortest path computations)
- src/LP_tools/preprocess.jl: clustering of data. Deprecated but kept because of potential usefulness.
- src/LP_tools/shortestpath.jl: copied from TaxiSimulation

#### Metropolis
- test/fast_LP_test.jl: wrapper function calls fast_LP on metropolis. fast_LP has changed a lot since this was last updated, so a fair amount of work needs to be done to make this work again. On the plus side, this will mostly require listing appropriate parameters
- test/loadMetropolis.jl: used to generate different metropolises, synthetic data
- test/metropolis_intro.jl: loads dependencies for metropolis

#### Visualization
- visualizer/manhattanVisualization.jl: code to draw Manhattan in SFML (helper)
- visualizer/algorithmOutput.jl: visualize algorithm results on Manhattan

## LP parameters
At the top of fast_LP.jl, there is a (very) long list of parameters to be fed into the method. This is inelegant but clear. Users should take time to familiarize themselves with these parameters. Parameters are specified in params.jl. Running params.jl will create a single parameter file that can then be used as input to the LP.

## Reading outputs

#### Output directories

Full Gurobi output logs are not saved, so output must be read from the src/Outputs/ directory.
The src/Outputs/ directory contains a directory for each test. Directory name specifies input parameters and is a little difficult to read at first. To see how the directory name is generated, take a look at the beginning of the fast_LP function in fast_LP.jl. Parameters for the run corresponding to an output directory are saved in parameters.json.

#### Contents of an output directory
- parameters.json: JSON file with parameters (same one used as input to qsub)
- info.txt: extremely outdated information file (deprecated by parameters.json). Only useful information is timestamp because it specifies when the test was run.
- timestats.csv: two-column CSV file, the first column is the iteration number, the second column is the runtime of that iteration in seconds. Last line only has a single number, which is the total runtime of all iterations.
- numODpairs.csv: two-column CSV file, the first column is the iteration number, and the second column is the number of origin-destination pairs in the model for that iteration.
- errorstats.csv: 7-column CSV file, these are what the columns stand for:

1. Iteration number
2. Total squared error (single-ride computation)
3. Average absolute relative error (single-ride computation)
4. Average bias (single-ride computation)
5. Total squared error (averaged node computation)
6. Average absolute relative error (averaged node computation)
7. Average bias (averaged node computation)
