# take a network as input
# simulate it 
# return output

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

# solveCSB(J)
# solveRacipe(J)

J = rand(-1:1, 3, 3)

D1 = calcFreq(solveCSB(J)) 
D2 = calcFreq(solveRacipe(J)) 

D1Si = parse.(Int, D2.Sequence)
bar(D1.Sequence, D1.RelFreq)
bar(D2.Sequence, D2.RelFreq)

function getMonoScores(J, evalfunc)
    D1 = calcFreq(evalfunc(J))
    S = getMonopositiveStrings(size(J, 1))
    f(x) = x in S 
    D1[f.(D1.Sequence), :]
    sort!(D1)
    D1 
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

    scatter(XCs', XRs')
end
=#