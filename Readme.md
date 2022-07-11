### Description

Collection of functions to start from an initial network (given or random) and use mutation and selection to arrive at a network with monopositive states, as evaluated with RACIPE. 

### Instructions

Clone the repository and change directory to `TopoFiles`. 

Run `julia --project`. 

Run `include("scriptEvolution.jl")`. 

Call function `simulateRacipe` with parameters in order:

1. A matrix (of integers) representing the network.
2. Number of iterations.
3. Number of nodes. _(only 3 supported now)_

Get a sequence of Matrices and their corresponding scores (which is the sum of relative frequencies of desired states). 
