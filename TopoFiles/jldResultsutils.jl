using CairoMakie

d = load("multiRacipeResults14-07-22-193136.jld")
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

"""Takes in the scoresmatrix `X` 
and plots the convergence over the iterations."""
function plotConvergence(X)
    plot(X', ylims=(0, 1),
        label="",
        color=:blue,
        lw=1, alpha=0.7,
        xlabel="iteration number",
        ylabel="score")
    xmean = mean(X, dims=1)
    plot!(xmean', label="mean score", color=:red, lw=2)
end

"""Takes in a matrix and returns its heatmap""" 
function plotNetworkHeatmap(network)
    ns = size(network,1)
    xs = 1:ns 
    ys = 1:ns
    CairoMakie.heatmap(xs, ys, reverse(network,dims=1))
    
end