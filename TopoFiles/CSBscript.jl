# take a network as input
# simulate it 
# return output
using Plots

"""constrains x between -1 and 1"""
function constrain(x, a::Float64=-1., b::Float64=1.)
    if x > b
        return b
    elseif x < a
        return a
    else 
        return x
    end
end

σ(x) = ((1/(1+exp(-10x)))-0.5)*2 

"""Take a network, evaluates it, returns the trajectory"""
function simulateCSBtraj(J)
    N = 100 #number of steps
    Δt = 0.1 #time step
    nn = size(J, 2) #number of nodes
    X = rand(-1:1.0:1, nn) #expression levels of each node 
    XN = zeros(nn, N)
    A = zeros(nn) 
    for k in 1:N
        A = J'X
        X = X + (Δt .* A) ./ (abs.(A) .+ 1)
        X = constrain.(X)
        XN[:, k] = X
    end 
    return XN
end

"""Testing the eval using sigmoid function"""
# function simulateCSBtest(J)
#     N = 100 #number of steps
#     Δt = 0.1 #time step
#     nn = size(J, 2) #number of nodes
#     X = rand(0:1.0:1, nn) #expression levels of each node 
#     XN = zeros(nn, N);
#     A = zeros(nn) 
#     for k in 1:N
#         A = J'X
#         X = X + (Δt .* A)
#         X = σ.(X)
#         XN[:, k] = X
#     end 
#     return XN
# end


"""Take a network, evaluates it, returns the trajectory"""
function simulateCSB(J)
    N = 1000 #number of steps
    Δt = 0.1 #time step
    nn = size(J, 2) #number of nodes
    X = rand(-1.:1., nn) #expression levels of each node 
    A = zeros(nn) 
    for _ in 1:N
        A = J'X
        X = X + (Δt .* A) ./ (abs.(A) .+ 1)
        X = constrain.(X)
    end 
    Xint = (x->constrain(x, 0., 1.)).(X) .|> round .|> Int .|> string
    return *(Xint...) 
end

solveCSB(J, n=1000) = [simulateCSB(J) for _ in 1:n]

function getMonoScores(J, evalfunc)
    D1 = calcFreq(evalfunc(J))
    S = getMonopositiveStrings(size(J, 1))
    f(x) = x in S 
    D1 = D1[f.(D1.Sequence), :]
    sort!(D1)
    D1 
end

function getAllBinaryStrings(n)
    x = String[]
    function f(s, x, n, i)
        if i < n
            f(s*"0", x, n, i+1)
            f(s*"1", x, n, i+1)
        else
            push!(x, s)
        end
    end
    return f("", x, n, 0)
end

#= Comparison between RACIPE and CSB =#

#=
begin
    Nt = 50
    XJ = Matrix[]
    XC = zeros(3, Nt)
    XR = zeros(3, Nt)

    for i in 1:Nt
        J = rand(-1:1, 3, 3)
        push!(XJ, J)
        XC[:, i] = getMonoScores(J, solveCSB).RelFreq
        XR[:, i] = getMonoScores(J, solveRacipe).RelFreq
    end

    heatmap(vcat(XC, XR),
        aspect_ratio=:equal,
        c=:viridis)

    XCs = sum(XC, dims=1)
    XRs = sum(XR, dims=1)

    #NORMALIZED
    XCn = (XCs .- mean(XCs)) ./ std(XCs)    
    XRn = (XRs .- mean(XRs)) ./ std(XRs)
    scatter(XCn', XRn', xlabel="CSB", ylabel="RACIPE", 
        legend=:none, 
        title="Sum of relative frequencies of monopositive states\nin RACIPE and CSB")
end
=#

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

function recombine(A::DataFrame, B::DataFrame, nnodes::Int)
    la, lb = nrow(A), nrow(B)
    #number of nodes
    getnn(A::DataFrame) = (cat(A.Source, A.Target, dims=1) |> unique |> length)
    nna = getnn(A) 
    nnb = getnn(B)
    if (nna != nnb)
        error("""The networks cannot be recombined, \
        not the same size (A: $(nna), B: $(nnb))""")
    end
    l = (la < lb) ? la : lb 
    nr = rand(1:l)
    Ae = A[nr+1:end, :]
    Be = B[nr+1:end, :]
    A = vcat(A[1:nr, :], Be)
    B = vcat(B[1:nr, :], Ae)
    removeduplicatepairs(A)
    removeduplicatepairs(B)
    nna = getnn(A) 
    nnb = getnn(B)
    if !(nna == nnb == nnodes)
        print("""the generated network size does not \
        match the specified size. Networks remain unmodified""") 
    end
    return A, B
end

nn = 3 #number of nodes
n = 6  #number of network pairs at each iteration
N = 5 #number of iterations 
X = [rand(-1:1, nn, nn) for i in 1:2n] #list of networks

# for i in 1:N
begin
    Xs = getPscore.(calcFreq.(solveCSB.(X)), "nrmsd") .+ 1
    X = sample(X, Weights(Xs), 2n, replace=true)
    for j in 1:2:2n
        X[j], X[j+1] = recombine(interaction2topo(X[j]), 
                                 interaction2topo(X[j+1]),
                                 nn) .|> topo2interaction
    end
end
