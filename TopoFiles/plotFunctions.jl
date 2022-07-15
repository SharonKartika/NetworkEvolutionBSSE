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



