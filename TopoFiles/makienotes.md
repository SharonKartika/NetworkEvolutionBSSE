### Figure

Main object. A canvas, to which we can add `Axis`, `Legend`, `Colorbar` etc.

`resolution`. Takes a tuple of width and height, respectively (in constrast to matrices).

Create axis in the figure, at position `1,1`   using `Axis(f[1,1])`. Can also add `title`, `xlabel`, `ylabel` by passing them as parameters.  (but don't do this. Gives weird plots).





```
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
```



```
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
```

