d = load("multiRacipeResults15-07-22-015055.jld")
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


