dist(A::Real, B::Real) = (A-B)^2 

function flatten(A)
    p = Any[]
    function flatten(A, p)
        if typeof(A) <: Real
            push!(p, A)
        else 
            for x in A
               flatten(x, p) 
            end
        end
    end
    flatten(A, p)
    return p
end

dist(A, B) = dist(mean(flatten(A)), mean(flatten(B)))

B = [1,2,3,7,8,9]

function hieclustvec(B)
    A = Any[]
    for i in B
        push!(A, i)
    end
    while length(A) > 2
        n = length(A)
        mindist = dist(A[1], A[2])  
        minargs = (1, 2)
        for i in 1:n-1
            for j in i+1:n
                d = dist(A[i], A[j])
                if d < mindist
                    mindist = d 
                    minargs = (i, j)
                end
            end
        end
        i, j = minargs 
        push!(A, [A[i], A[j]])
        deleteat!(A, minargs)
    end
    return A
end