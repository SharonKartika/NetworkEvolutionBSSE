include("../Boolean.jl/bmodel.jl")
include("../Boolean.jl/utils.jl")
using StatsBase
using DataFrames, CSV, Statistics
using Plots
using JLD
using GraphRecipes
using Suppressor
using Dates


function mutate(network::Array{Int64,2}, score::Float64, nr::Int)
    """Accepts a network (as a matrix), its fitness score (in [0,1]), and number of elements to replace `nr`;
    Returns a mutated network depending on the score;
    High fitness: low mutation; Low fitness: high mutation"""

    # replaces `nr` elements in `network` randomly 
    network[sample(1:9, nr, replace=false)] = rand(-1:1, nr)
    return network
end

function mutateMulti(network::Array{Int64,2}, score::Float64, n::Int)
    """Takes a network and a fitness score
    Returns `n` mutated matrices"""

    #number of elements in the network
    nm = size(network)[1]
    nm2 = nm^2
    # maps score (in [0, 1]) to 1:9
    # nr = Int(round(nm2 -score*(nm2-1))) #nr depends on score 
    nr = 3 # hardcoded nr. nr elements are mutated always

    # vector of `n` matrices of size `nm x nm`
    # find a faster way
    nNets = [copy(network) for i in 1:n]
    # vectorized mutation
    mutate.(nNets, score, nr)
end

function getPscore(df)
    "returns the total probability of getting any of the states in `vars`; takes in the frequency df"
    vars = ["'001'"; "'010'"; "'100'"]
    sum(df[ainb(df.states, vars), :frequency])
    # freqList = df[ainb(df.states,vars),:frequency] 
    # if length(freqList) < 3  
    #     return 0.
    # else
    #     return prod(freqList)
    # end
end

function dfFreqMulti(state_df::Array{DataFrame,1}, cols::Array{Symbol,1})
    "returns an array of dataframes containing the frequencies of states in each network "
    return [dfFreq(i, [:fin, :flag]) for i in state_df]
end

function ainb(a, b)
    "find index of a in b"
    [i for i in 1:length(a) if a[i] in b]
end

function interaction2topo(tNet::AbstractMatrix, fname::Int)
    xi, yi, vi = findnz(sparse(tNet))
    lets = 'A':'Z'
    nmap(x) = (3 - x) ÷ 2 # to map 1->1 and -1->2
    df = DataFrame(Source=lets[xi], Target=lets[yi], Type=nmap.(vi))
    CSV.write("net$(fname).topo", df, delim='\t')
    df
end

toggleTriad = [0 -1 -1; -1 0 -1; -1 -1 0]

network2 = [1 -1 -1; -1 0 -1; -1 -1 1];



function solveRacipe(network::Matrix{Int64}, i::Int64)
    "Takes a network, and a reference number as input. Returns the frequency matrix of states"
    interaction2topo(network, i)
    output = @capture_out run(`./RACIPE net$(i).topo -threads 40`)
    dfr = CSV.read("net$(i)_solution.dat", DataFrame; header=0)
    n = (dfr |> names |> length) ÷ 2
    for i = 1:n
        meani = mean(dfr[:, n+i])
        dfr[:, "Column$(2n+i)"] = (dfr[:, n+i] .- meani) .> 0.0
    end
    dfr[:, "fin"] = "'" .* (.*([string.(Int.(dfr[:, 2n+i])) for i in 1:n]...)) .* "'"
    dfFreq(dfr, [:fin])
end

function simulateBoolean(niter=10, nnodes=3)
    #step0: storing scores and matrices 
    x = zeros(niter)
    Mi = Array{Float64}(undef,
        nnodes,
        nnodes,
        niter)#best network in each iteration 
    println("0%|$("-"^niter)|100%") #progress
    print("   ")

    # step1: loading and initial evaluation 
    initNetwork = rand((-1:1), (nnodes, nnodes))
    df = dfFreq(
        asyncUpdateStates(initNetwork, 10000, 1000),
        [:fin, :flag])
    pscore = getPscore(df)

    for i in 1:niter
        # step2: mutate
        nNetworks = mutateMulti(initNetwork, pscore, 6)
        push!(nNetworks, initNetwork)

        # step3: evaluate 
        nResults = asyncUpdateStates.(nNetworks, 1000, 1000) # find states of all networks
        nFreqs = dfFreqMulti(nResults, [:fin, :flag]) # find frequencies of states (VECTORIZE)

        # step4: select
        pscores = getPscore.(nFreqs) # find pscores of all networks
        topScoreIndex = argmax(pscores) # find index of highest scoring network 
        initNetwork = nNetworks[topScoreIndex] # select the highest scoring network

        # display(pscores[topScoreIndex]), display(initNetwork)
        x[i] = pscores[topScoreIndex]
        Mi[:, :, i] = copy(initNetwork)
        print("=")
    end
    print("\n")
    return x, Mi
end


function simulateRacipe(initNetwork::Matrix{Int64}, niter=10, nnodes=3)
    #step0: storing scores and matrices 
    x = zeros(niter)
    Mi = Array{Int64}(undef,
        nnodes,
        nnodes,
        niter)#best network in each iteration 
    println("0%|$("-"^niter)|100%") #progress
    print("  |")
    # step1: loading and initial evaluation 
    # initNetwork = rand((-1:1), (nnodes, nnodes));
    df = solveRacipe(initNetwork, 1)
    pscore = getPscore(df)

    timeTaken = @elapsed for i in 1:niter
        # step2: mutate
        nNetworks = mutateMulti(initNetwork, pscore, 6)
        push!(nNetworks, initNetwork)

        # step3: evaluate 
        # find states of all networks 
        nFreqs = [solveRacipe(net, 1) for net in nNetworks]

        # step4: select
        pscores = getPscore.(nFreqs) # find pscores of all networks
        topScoreIndex = argmax(pscores) # find index of highest scoring network 

        initNetwork = nNetworks[topScoreIndex] # select the highest scoring network
        # display(pscores[topScoreIndex]), display(initNetwork)
        x[i] = pscores[topScoreIndex]
        Mi[:, :, i] = copy(initNetwork)
        # println("$(i)%")
        print("=")
    end
    print("|\n$(timeTaken) seconds\n")
    return x, Mi
end


#=begin # testing pscores
    df = DataFrame(states=["'001'", "'011'", "'010'", "'100'"],
                frequency=[1.,0. ,2., 3.])
    df = DataFrame(states=["'001'", "'011'"],
            frequency=[4.2,2.])

    getPscore(df)
end=#


#=begin #test boolean toggleTriad
    x=Float64[]
    for i in 1:300
    df = dfFreq(
            asyncUpdateStates(toggleTriad, 1000, 1000),
            [:fin, :flag])
        pscore = getPscore(df)
        push!(x, pscore)
    end
    plot(x, ylims=[0., 0.007])
end=#

#old interaction2topo. Relies on manual actions
#=function interaction2topo(tNet::AbstractMatrix, fname::Int)
    xi, yi, vi = findnz(sparse(tNet))
    lets = 'A':'Z'
    nmap = Dict(1=>1, -1=>2)
    io = open("net$(fname).topo", "w") 
    write(io,"Source\tTarget\tType\n")
    for i in eachindex(vi)
        # println("$(lets[xi[i]])\t$(lets[yi[i]])\t$(nmap[vi[i]])")  
        write(io, "$(lets[xi[i]])\t$(lets[yi[i]])\t$(nmap[vi[i]])\n")
    end
    close(io)
end=#

#= #an ensemble of trajectories  
# will take 50*20*7*2.5/60 = 291 minutes
X = Any[]
M = Any[]
for i in 1:20
    x, m = simulateRacipe(50, 3)
    push!(X, x)
    push!(M, m)
end
=#

# start with toggleTriad and simulate 

# x, Mi = simulateRacipe(toggleTriad, 5, 3)



function multiRacipe(network)
    # start with the network passed for all iterations
    niter = 100 # number of iterations of 
    nrepl = 20 # number of replicates
    scoresMatrix = Array{Float64,2}(undef, nrepl, niter)
    curdate = Dates.format(Dates.now(), "dd-mm-yy-HHMMSS")
    networkList = []
    for i in 1:nrepl
        print("$(i)/$(nrepl)\n")
        x, Mi = simulateRacipe(network, niter, 3)
        scoresMatrix[i, :] = x
        push!(networkList, Mi)
    end
    save("multiRacipeResults$(curdate).jld", "scoresMatrix", scoresMatrix, "networkList", networkList)
    return scoresMatrix, networkList
end

function multiRacipe()
    # start with random network for each iteration if no network is passed 
    niter = 40 # number of iterations of 
    nrepl = 15 # number of replicates
    scoresMatrix = Array{Float64,2}(undef, nrepl, niter)
    curdate = Dates.format(Dates.now(), "dd-mm-yy-HHMMSS")
    networkList = []
    for i in 1:nrepl
        print("$(i)/$(nrepl)\n")
        network = rand(-1:1, 3, 3)
        x, Mi = simulateRacipe(network, niter, 3)
        scoresMatrix[i, :] = x
        push!(networkList, Mi)
    end
    save("multiRacipeResultsRandomInitial$(curdate).jld", "scoresMatrix", scoresMatrix, "networkList", networkList)
    return scoresMatrix, networkList
end

function plotMultiRacipeConvergence(network)
    X, XM = multiRacipe(network)
    plot(X', ylims=[0, 1], legend=:none)
end


