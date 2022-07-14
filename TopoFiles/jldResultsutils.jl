d = load("multiRacipeResults14-07-22-193136.jld")
X, XM = d["scoresMatrix"], d["networkMatrix"]

function getNetworkScore(nthrun, nthiter, X, XM)
    i = nthrun
    j = nthiter
    XM[i][:, :, j], X[i, j]
end

function getNetworkScore(indice, X, XM)
    getNetworkScore(indice[1], indice[2], X, XM)
end

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

    selectedNetworks = []
    for i in 1:noverthresh
        p = getNetworkScore(indices[i, :], X, XM)[1]
        if !(p in selectedNetworks)
            push!(selectedNetworks, p)
        end
    end
    selectedNetworks
end

function plotMatrixHeatmap(network)
    heatmap(reverse(network, dims=1),
        colorbar=:none,
        aspect_ratio=1,
        showaxis=false,
        ticks=false,
        color=[:turquoise3, :wheat1, :tomato2])
end

function matrixSubplots(selectedNetworks)
    l = @layout grid(5, 8)
    plot([plotMatrixHeatmap(selectedNetworks[i]) for i in 1:length(selectedNetworks)]..., layout=l)
end

function plotConvergence(X)
    plot(X', ylims=(0,1),
    label="",
    color=:blue,
    lw=1, alpha=0.7,
    xlabel="iteration number",
    ylabel="score")
    xmean = mean(X, dims=1)
    plot!(xmean', label="mean score", color=:red, lw=2)
end
