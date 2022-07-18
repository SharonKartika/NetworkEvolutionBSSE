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

        #directions correct; colors wrong
        add_edge!(g, dfe[i, "SourceInt"], dfe[i, "TargetInt"])
        
        # colors correct, directions wrong
        # add_edge!(g, dfe[i, "TargetInt"], dfe[i, "SourceInt"])
        
    end
    
    edgeColors = mapcolors.(dfe.Type)
    f, ax, p = GraphMakie.graphplot(g;
        arrow_size=25,
        arrow_shift=0.5,
        edge_color=edgeColors,
        nlabels=string.(intToLetter.(1:size(network, 1))))

    hidedecorations!(ax); hidespines!(ax); ax.aspect = DataAspect()
    return f
end

"""Takes in a vector of matrices, and plots their graphs"""
function plotNetworkMulti(networks)
    letterToInt(letter) = (Dict('A':'Z' .=> 1:26))[letter]
    intToLetter(i) = (Dict(1:26 .=> 'A':'Z'))[i]
    mapcolors(i) = ((i==2) ? (:red) : (:green))

    fig = Figure()
    count = 1
    axs = []
    for i in 1:2
        for j in 1:5
            network = networks[count]
            g = SimpleDiGraph(size(network,1))

            dfe = interaction2topo(network)
            dfe[:, "SourceInt"] = letterToInt.(dfe[:,"Source"])
            dfe[:, "TargetInt"] = letterToInt.(dfe[:, "Target"]) 
        
            for k in 1:nrow(dfe)
                add_edge!(g, dfe[k, "SourceInt"], dfe[k, "TargetInt"])
            end
            edgeColors = mapcolors.(dfe.Type)

            ax = Axis(fig[i,j], aspect=DataAspect())
            push!(axs, ax)
            GraphMakie.graphplot!(ax, g;
                arrow_size=20,
                arrow_shift=0.9,
                edge_color=edgeColors,
                nlabels=string.(intToLetter.(1:size(network, 1))))
    
            count+=1
        end
    end
    hidedecorations!.(axs); hidespines!.(axs);
    fig
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

  """Takes in the scoresmatrix `X` 
  and plots the convergence over the iterations."""
  function plotConvergenceTest(X)
      fig, ax, sp = series(X,
                      solid_color=(:blue,0.4))
                      # labels=[:nothing for i in 1:size(X,1)])
      ylims!(ax, (0, 1))
      xlims!(ax, (0, size(X,2)))
      lp = lines!(ax, vec(mean(X,dims=1)), color=:red,linewidth=3,label="mean score")
      ax.xlabel = "Iterations"
      ax.ylabel = "Score"
      ax.title = "Convergence curve for $(size(X,1)) runs\nand $(size(X,2)) iterations"
      axislegend(ax,
                  [sp, lp],
                  ["score for each run", "mean score"] )
      fig
  end
