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
N = influenceMatrix(N)
# heatmap(N, aspect_ratio=:equal,
#         c=:viridis)

# M = influenceMatrix(N)

dist(A, B) = sum((A .- B) .^ 2)
nn = size(N, 1)

Dc = zeros(nn, nn)
Dr = zeros(nn, nn)
for i in 1:nn-1
    for j in i+1:nn
        Dc[i, j] = Dc[j, i] = dist(N[:, i], N[:, j])
        Dr[i, j] = Dr[j, i] = dist(N[i, :], N[j, :])        
    end
end

hccol = hclust(Dc)
hcrow = hclust(Dr)

ctc2 = cutree(hccol, k=2)
ihc1 = [i for i in 1:length(ctc2) if ctc2[i]==1]
ihc2 = [i for i in 1:length(ctc2) if ctc2[i]==2] 
ctr2 = cutree(hcrow, k=2)
ihr1 = [i for i in 1:length(ctr2) if ctr2[i]==1]
ihr2 = [i for i in 1:length(ctr2) if ctr2[i]==2] 

M = copy(N)
M = [M[:, ihc1]  M[:, ihc2]]
M = [M[ihr1, :]; M[ihr2, :]]

heatmap(M, aspect_ratio=:equal, cbar=:none, c=:RdBu)
heatmap(M, aspect_ratio=:equal, cbar=:none, c=:thermal)

