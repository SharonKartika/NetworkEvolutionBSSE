using CairoMakie,
    GraphMakie,
    Graphs,
    GraphMakie.NetworkLayout

include("scriptEvolutionGeneral.jl")

function getGraph()

    g = SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 4)
    add_edge!(g, 4, 3)
    add_edge!(g, 3, 2)
    add_edge!(g, 2, 5)
    add_edge!(g, 5, 4)
    add_edge!(g, 4, 1)
    add_edge!(g, 1, 5)

    # define some edge colors
    edgecolors = [:black for i in 1:ne(g)]
    edgecolors[4] = edgecolors[7] = :red

    f, ax, p = GraphMakie.graphplot(g, layout=Shell(),
        node_color=[:black, :red, :red, :red, :black],
        edge_color=edgecolors)

    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()
    return f
end

function graphPlotSelf()
    g = complete_graph(3)
    add_edge!(g, 1, 1)
    add_edge!(g, 2, 2)
    add_edge!(g, 3, 3)
    f, ax, p = GraphMakie.graphplot(g)

    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()
    return f
end

"""Takes in a interaction matrix,
 and plots the network corresponding to it"""
function plotNetwork(network)
    letterToInt(letter) = (Dict('A':'Z' .=> 1:26))[letter]
    mapcolors(i) = ((i==2) ? (:red) : (:green))
    g = SimpleDiGraph(size(network,1))
    
    dfe = interaction2topo(network)
    dfe[:, "SourceInt"] = letterToInt.(dfe[:,"Source"])
    dfe[:, "TargetInt"] = letterToInt.(dfe[:, "Target"]) 
    
    for i in 1:nrow(dfe)
        add_edge!(g, dfe[i, "SourceInt"], dfe[i, "TargetInt"])
    end
    
    edgeColors = mapcolors.(dfe.Type)
    f, ax, p = GraphMakie.graphplot(g; arrow_size=0,
        edge_color=edgeColors)
    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()
    return f
end



