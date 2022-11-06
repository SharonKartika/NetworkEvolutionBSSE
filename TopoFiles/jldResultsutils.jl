d = load("multiRacipeResults15-07-22-015055Random4by4.jld")
X, XM = d["scoresMatrix"], d["networkMatrix"]

"""Returns the score and the matrix of the network
 in `nthiter` of `nthrun`."""
function getNetworkScore(nthrun, nthiter, X, XM)
    X[nthrun, nthiter], XM[nthrun, nthiter]
end

"""Returns the score and the matrix of the network
 in `nthiter` of `nthrun`. Takes in a tuple of indices"""
function getNetworkScore(indice, X, XM)
    getNetworkScore(indice[1], indice[2], X, XM)
end

"""Returns a vector of networks which have scores
 above the threshold `thresh`"""
function findNetworksAboveThresh(thresh, X, XM)
    noverthresh = sum(X .> thresh)
    indices = zeros(Int, noverthresh, 2)

    count = 1
    for i in 1:size(X)[1], j in 1:size(X)[2]
        if X[i, j] > thresh
            indices[count, :] = [i, j]
            count += 1
        end
    end

    selectedNetworks = Vector{Matrix{Int}}(undef, 0)
    for i in 1:noverthresh
        p = getNetworkScore(indices[i, :], X, XM)[2]
        if !(p in selectedNetworks)
            push!(selectedNetworks, p)
        end
    end
    selectedNetworks
end

"""Solves a network, 
finds the frequency of each monoposotive state"""
function getDictFreq(network)
    sols = solveRacipe(network)
    # D = countmap(sols)
    D = proportionmap(sols) .= 0
    S = getMonopositiveStrings(size(network, 1))
    Dnew = Dict()
    for i in S
        if i in keys(D)
            Dnew[i] = D[i]
        end
    end
    Dnew
end

"""Barplots the relative frequency of monopositive states.
Takes the output of getDictFreq"""
function barplotDictFreq(dfreq)
    x = String[]
    y = Float64[]
    for i in keys(dfreq)
        push!(x, i)
    end
    for i in values(dfreq)
        push!(y, i)
    end

    fig = barplot(y, axis=(; xticks=(1:length(x), x)))
    ylims!(0, 1)
    fig
end


function getMonoFreq(network)
    """returns the frequencies of monopositive strings
     for a network"""
    D = calcFreq(solveRacipe(network))
    nn = length(D[1, 1])
    S = getMonopositiveStrings(nn)
    Rfs = filter(:Sequence => x -> x in S, D)
    sort!(Rfs)
    return Rfs.RelFreq
end

function plotMonoPositiveMulti(D,
     xtext::String="net",
     plottitle::String="""Heatmap of \
     relative frequencies for different networks\n  """)
    
    "Takes the output of `getMonoFreq`
     applied on a set of networks, and plots
     the frequencies of monopositive states,
     as a heatmap"
    # D = getMonoFreq.(networks)
    Dm = reduce(hcat, D) #conv. vector{vector} to Matrix
    nn = size(Dm,1)
    xx = size(Dm, 2)
    fig = Figure()
    ax = Axis(fig[1,1],
         aspect=DataAspect(),
         xticks=(1:xx, (x->"$(xtext) $x").(1:xx)),
         yticks=(1:nn,getMonopositiveStrings(nn))
         )
    x = size(Dm, 2)
    hm= heatmap!(ax, Dm',
         colorrange=(0,1/nn))
    Colorbar(fig[2,:], hm, vertical=false)
    ax.title=plottitle
    fig
end


