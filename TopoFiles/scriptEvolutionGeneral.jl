#=
Genetic algorithms to select networks (of arbitrary sizes)
giving monopositive states size using RACIPE
=#

using StatsBase,
    DataFrames,
    CSV,
    Statistics,
    Plots,
    JLD,
    GraphRecipes,
    Suppressor, 
    Dates, 
    SparseArrays,
    LinearAlgebra

toggleTriad = [0 -1 -1; -1 0 -1; -1 -1 0]
toggleSquare = (-1)*(ones(Int, 4,4)- I)

"""Takes a network, returns a mutant with `nr` elements replaced"""
function mutate(network, nr)
    network[sample(1:length(network), nr, replace=false)] = rand(-1:1, nr)
    return network
end

"""Takes a network, returns `n` mutated networks,
including the unmodified original network"""
function mutateMulti(network, n)
    nNets = [copy(network) for i in 1:n]
    mutate.(nNets[2:end], 5) 
    return nNets
end

"find indices of elements of `a` that are present in `b`"
function ainb(a, b)
    [i for i in 1:length(a) if a[i] in b]
end

"""Returns a list of strings containing `1`s and `0`s,
which are monopositive"""
function getMonopositiveStrings(n)
    [("0"^(i-1))*"1"*("0"^(n-i)) for i in 1:n]
end

"""Takes in the df of frequencies of occurrence (`dfFreq`), 
and calculates the score of the network"""
function getPscore(dfFreq)
    nn = length(dfFreq[1,1]) #number of nodes in network
    reqStates = getMonopositiveStrings(nn) # list of required states
    indexOfReq = ainb(dfFreq.Sequence, reqStates)
    return sum(dfFreq[indexOfReq, "RelFreq"])
end

"""Takes in the result of simulation
(a list of strings of outputs), and finds the relative
frequencies of each unique state."""
function calcFreq(dfr)
    D = proportionmap(dfr)
    dfFreq = DataFrame(Sequence=collect(keys(D)), RelFreq=collect(values(D)))
end

"""Takes an interaction matrix, converts it to topo format,
and returns it"""
function interaction2topo(tNet::AbstractMatrix)
    xi, yi, vi = findnz(sparse(tNet))
    lets = 'A':'Z'
    nmap(x) = (3 - x) ÷ 2 # to map 1->1 and -1->2
    df = DataFrame(Source=lets[xi], Target=lets[yi], Type=nmap.(vi))
    CSV.write("netevol.topo", df, delim='\t')
    df
end

"""Runs RACIPE on the network. Returns a vector of final states."""
function solveRacipe(network)
    interaction2topo(network)

    #run RACIPE; store RACIPE output in output
    output = @capture_out run(`./RACIPE netevol.topo -threads $(Threads.nthreads())`)
    
    #read the results; ignore first 3 columns
    dfr = CSV.read("netevol_solution.dat", DataFrame; header=0)[:, 4:end]
    n = dfr |> names |> length 
    
    rename!(dfr, ["Column$(i)" for i in 1:n])
    for i in 1:n
        meani = mean(dfr[:, i])
        #if value greater than mean of the column, set to 1; else 0
        dfr[:, "Column$(n+i)"] = string.(Int.(dfr[:, i] .> meani))
    end
    
    #join the state of individual nodes and store as string
    dfr[:, "fin"] = .*([dfr[:, i] for i in n+1:2n]...)
    dfr[:, "fin"]
end

"""Runs the evolutionary algorithm for `niter` iterations
Returns p-score as well as the best network at each iteration.    """
function simulateRacipe(network, niter=10)
    #step0: storing scores and matrices 
    nnodes = size(network, 1)
    x = zeros(niter) #score of best network
    # best network in each iteration
    Mi = Vector{Matrix{Int64}}(undef, niter)

    #display progress
    println("0%|$("-"^niter)|100%") 
    print("  |")

    # step1: initial evaluation 
    #pscore = solveRacipe(network) |> calcFreq |> getPscore

    timeTaken = @elapsed for i in 1:niter
        # step2: mutate
        nNetworks = mutateMulti(network, 6)
        # step3: evaluate 
        # find states of all networks 
        nResults = [solveRacipe(net) for net in nNetworks]
        nDfFreqs = calcFreq.(nResults)
        pscores = getPscore.(nDfFreqs) # find pscores of all networks

        # display((nNetworks .=> pscores))

        # step4: select and store
        topScoreIndex = argmax(pscores) # find index of highest scoring network 
        network = nNetworks[topScoreIndex] # select the highest scoring network

        x[i] = pscores[topScoreIndex]
        Mi[i] = copy(network)

        print("=")
    end # step5: go to step2 if i < niter
    print("|    $(timeTaken) seconds\n")
    return x, Mi
end

