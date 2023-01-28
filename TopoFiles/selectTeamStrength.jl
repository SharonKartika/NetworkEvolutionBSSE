J = rand(-1:1, 10, 10)


"""Scoring function. Takes """
function costFunction(J::AbstractMatrix)
    return sum(J)    
end

function best(M, costFunction::Function)
    scores = costFunction.(M)
    i = argmax(scores)
    return M[i]
end

function createMutant(J::AbstractMatrix; nr::Int = 0)
    l = length(J)
    (nr==0) ? (nr = l ÷ 10) : () #10% mutation rate 
    J[sample(1:l, nr, replace=false)] = rand(-1:1, nr)
    # ensure that all the nodes are connected.  
    # see if the network has empty crosses. If yes, add random element in the row.
    # nzeros = zeros(Int, n)
    # for i in 1:n
    #     if (J[i, :] == nzeros) && (J[:, i] == nzeros)
    #         J[i, rand(1:n)] = rand([-1, 1])
    #     end
    # end
    return J
end

function createMutants(J::AbstractMatrix;
                       nmutants::Int = 6,
                       nr::Int = 0)
    return createMutant.([copy(J) for _ in 1:nmutants], nr=nr)
end

function evolveNet(J::AbstractMatrix,
                   costFunction::Function;
                   nmutants::Int=10,
                   nr::Int = 0,
                   niter::Int = 100)
    C = copy(J)
    seen = [C]
    selected = [C]
    for i in 1:niter
        if (i%50 == 0)
            println(i)
        end
        M = createMutants(C, nmutants=nmutants, nr=nr)
        M = setdiff(M, seen)
        if length(M) == 0
            continue
        end
        C = best(M, costFunction)
        push!(seen, M...)
        push!(selected, C)
    end
    return selected
end

