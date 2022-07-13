### Evolving Gene Regulatory Networks

We simulate gene regulatory networks using RACIPE and use a genetic algorithm to select GRNs that produce desired outcomes. 

### Instructions

Clone the repo. Change directory to `TopoFiles`.

Run julia (`julia --project`).

Include the required program using `include("scriptEvolutionGeneral.jl")`

Call `simulateRacipe` passing as arguments the network and the number of iterations respectively. The network is represented by a square matrix of integers -1, 0, or 1. The function returns the best network and its p-score, obtained at each iteration. The p-score is the sum of relative frequencies of the desired outcomes.  		

 
