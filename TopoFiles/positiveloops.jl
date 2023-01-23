using Graphs

nn = 3
N = rand(-1:1, nn, nn)

"""Takes an interaction matrix, converts it to a list of edges
Loses activation/inhibition information"""
function interaction2edgelist(tNet::AbstractMatrix)
    xi, yi, vi = findnz(sparse(tNet))    
    return [(xi[i], yi[i]) for i ∈ eachindex(xi)], vi
end

E, F = interaction2edgelist(N)#effect
E = Edge.(E)
G = SimpleDiGraph(E)

C = simplecycles_limited_length(G, 10)

#find number of positive loops

pfbl = 0
nfbl = 0
for loop in C
    v = 1
    for i in loop
        v *= F[i]
    end
    if v == 1
        pfbl += 1
    else 
        nfbl += 1
    end
end
pfbl, nfbl