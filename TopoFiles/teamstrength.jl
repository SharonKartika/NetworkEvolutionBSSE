using Clustering
using Plots 

function influenceMatrix(J::Matrix, lmax::Int = 10)
    nonzerodiv(x, y) = (y==zero(y)) ? (zero(y)) : (x / y) 
    Adjˡ = copy(J)
    A = Float64.(Adjˡ .!= 0)
    Adjˡmax = copy(A)
    Infl_max = zero(A)
    for l in 2:lmax
        Adjˡmax = Adjˡmax * A 
        Adjˡ = Adjˡ * J     
        Infl_max .+= (nonzerodiv.(Adjˡ, Adjˡmax))   
    end
    Infl_max ./= lmax
    return Infl_max
end

N = rand(-1:1, 10, 10)
heatmap(N, aspect_ratio=:equal,
        c=:viridis)
    
M = hclust(N+N',linkage=:average)