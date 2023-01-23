# Evolving gene regulatory networks



### Instructions

Clone the repository. Change directory into `TopoFiles`. 

Run julia as project `julia --project`. Optionally, specify number of threads.

Include the necessary files using `include("scriptEvolutionGeneral.jl")`. 

To run the algorithm, call `simulateRacipe` with the required arguments (see docstring in the code file). *Remeber to capture the return values*.

To do multiple runs from the same initial network, call `multiRacipe` with the required arguments.  Automatically stores the results in a file with `.jld` extension (see extras).



### Extras

#### Visualising networks

Include the necessary files using `include"plotFunctions.jl"`.

Call `plotNetwork` to plot the GRN as a graph. See docstring in code file for arguments.

Call `plotNetworkHeatmap` to plot the GRN as a heatmap of the interaction matrix.

Multiple networks can be plotted in the same plot using the functions `plotNetworkHeatmapMulti` and `plotNetworkMulti`. See docstrings in the code file.

To plot the convergence pattern of the result of `multiRacipe`, call `plotConvergence` with the required arguments.



#### Reading `.jld` files

Open the `jldResultsutils.jld` file and edit the first line to the required file.

Include the necessary files using `include("jldResultsutils.jld")`.

The results of the `multiRacipe` run is now stored in variables `X` (the matrix of scores), and `XM` (the matrix of networks). 



