include("scriptEvolutionGeneral.jl")

n = 6 #number of candidates
nn = 3 #number of nodes 
t = 4 #number of iterations 
startingNetworks = [rand(-1:1, nn, nn) for i in 1:n]
candidateNetworks = copy(startingNetworks)
candidateScores = zeros(n)

#=
TODO 
- [x] Filter out nn < 3
- [x] nRMSD lies between 0, -1. Map it to 0, 1;
- [] add mutations
- [] Check that recombine works
=#

# for t in 1:t
# populate in proportion to score 
candidateScores = getPscore.(calcFreq.(solveRacipe.(candidateNetworks)), "nrmsd");
candidateScores .+= 1 #map from (-1, 0) to (0, 1);
candidateNetworks = sample(candidateNetworks,
                            Weights(candidateScores),
                            n, 
                            replace=true)
# recombine 
candidateTopofiles = interaction2topo.(candidateNetworks)

"""Eliminates duplications. Modifies original dataframe"""
function removeduplicatepairs(A::DataFrame)
    n = nrow(A)
    z = [(A[i, :Source], A[i, :Target]) for i in 1:n] #array of pairs
    duplicateindex = Int[]
    for i in 1:length(z)-1
        for j in i+1:length(z)
            if z[i]==z[j]
                push!(duplicateindex, j)
            end
        end
    end
    deleteat!(A, duplicateindex)
end

function recombine(
                    A::DataFrame,
                    B::DataFrame,
                    nnodes::Int #number of nodes in matrix repr.
                  )
    #number of genes (rows in topofile)
    la, lb = nrow(A), nrow(B)
    #number of nodes
    getnn(A::DataFrame) = (cat(A.Source, A.Target, dims=1) |> unique |> length)
    nna = getnn(A) 
    nnb = getnn(B)
    if (nna != nnb)
        error("""The networks cannot be recombined, \
        not the same size (A: $(nna), B: $(nnb))""")
    end
    #shorter length 
    l = (la < lb) ? la : lb 
    #position to recombine
    nr = rand(1:l)
    #rows after the position to recombine. Potentially empty  
    Ae = A[nr+1:end, :]
    Be = B[nr+1:end, :]
    #recombination
    A = vcat(A[1:nr, :], Be)
    B = vcat(B[1:nr, :], Ae)
    ## size checks
    # uniqueness
    removeduplicatepairs(A)
    removeduplicatepairs(B)
    # number of nodes 
    nna = getnn(A) 
    nnb = getnn(B)
    if !(nna == nnb == nnodes)
       print("""the generated network size does not \
       match the specified size. Networks remain unmodified""") 
    end
    # return A, B 
end

