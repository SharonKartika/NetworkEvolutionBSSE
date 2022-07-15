using CairoMakie,
    GraphMakie,
    Graphs,
    GraphMakie.NetworkLayout

include("scriptEvolutionGeneral.jl")


"""Takes in a interaction matrix,
 and plots the network corresponding to it"""
function plotNetwork(network)
    letterToInt(letter) = (Dict('A':'Z' .=> 1:26))[letter]
    intToLetter(i) = (Dict(1:26 .=> 'A':'Z'))[i]
    mapcolors(i) = ((i==2) ? (:red) : (:green))
    g = SimpleDiGraph(size(network,1))
    
    dfe = interaction2topo(network)
    dfe[:, "SourceInt"] = letterToInt.(dfe[:,"Source"])
    dfe[:, "TargetInt"] = letterToInt.(dfe[:, "Target"]) 
    
    for i in 1:nrow(dfe)
        add_edge!(g, dfe[i, "SourceInt"], dfe[i, "TargetInt"])
    end
    
    edgeColors = mapcolors.(dfe.Type)
    f, ax, p = GraphMakie.graphplot(g;
        arrow_size=20,
        arrow_shift=0.9,
        edge_color=edgeColors,
        nlabels=string.(intToLetter.(1:size(network, 1))))

    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()
    return f
end

"""Takes a matrix and plots its heatmap"""
function plotNetworkHeatmap(network)
    fig = Figure()
    ax = Axis(fig[1,1], aspect=DataAspect())
    revnet = reverse(network,dims=2)
    nsize = size(network,2)
    heatmap!(ax, revnet)
    text!(ax,
         [string(revnet[i,j]) for i in 1:nsize for j in 1:nsize ],
         position=[(x, y) for x in 1:nsize for y in 1:nsize],
         textsize=30,
         color=[revnet[i,j]<0 ? (:white) : (:black) for i in 1:nsize for j in 1:nsize],
         align=(:center, :center))

    hidedecorations!(ax)
    hidespines!(ax)
    fig
end

"""Takes a vector of matrices and plots their heatmaps"""
function plotNetworkHeatmapMulti(networks)
    fig = Figure()
    count = 1
    axs = []
    for i in 1:2
        for j in 1:5
            revnet = reverse(networks[count],dims=2)
            nsize = size(revnet,2)
            ax = Axis(fig[i,j], aspect=DataAspect())
            push!(axs, ax)
            heatmap!(ax, revnet)
            text!(ax,
              [string(revnet[i,j]) for i in 1:nsize for j in 1:nsize ],
              position=[(x, y) for x in 1:nsize for y in 1:nsize],
              textsize=15,
              color=[revnet[i,j]<0 ? (:white) : (:black) for i in 1:nsize for j in 1:nsize],
              align=(:center, :center))
            count+=1
        end
    end
    hidedecorations!.(axs)
    hidespines!.(axs)
    fig
  end

