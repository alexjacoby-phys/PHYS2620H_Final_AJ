using LinearAlgebra
using Combinatorics
using Plots


#creates list of operators to pass to makeop#
function MakeOpList(A::Array{Complex{Float64},2},B::Array{Complex{Float64},2},x::Int64, y::Int64,L::Int64)
    Id = [1. 0 0 ; 0 1 0; 0 0 1]
    OpList = Vector{Array{Complex{Float64},2}}(undef,L)
    for i in 1:L
        if i == x
            OpList[i]= A
        elseif i == y
            OpList[i] = B
        else
            OpList[i] = Id
        end
    end
    return OpList
end


function makerho(statein::Array{Complex{Float64},2},n::Int64)
    rhostate = statein[n,:]
    rho = Array{Complex{Float64},2}
    rho = kron(rhostate,transpose(rhostate))
    rho = (1/tr(rho))*rho
    return rho
end





#does tensor product to make operator act on the whole space#
function MakeOp(C::Vector{Array{Complex{Float64},2}})
    HTerm = C[1]
    for i in 2:length(C)
        HTerm = kron(HTerm,C[i])
    end
    return HTerm
end

function traceop(n)
    trop = [1. 0 0 ; 0 1 0 ; 0 0 1]*(1+0*im)
    for i in 1:n-2
        trop = kron(trop,id)
    end
    trop = kron(trop,[1 1 1])
    return trop
end

function partialtrace(inputrho::Array{Complex{Float64},2},k)
    for i in 1:k
        inputrho = traceop(L+1-i)*inputrho*transpose(traceop(L+1-i))
    end
    return inputrho
end

function SvN(dmat::Array{Complex{Float64},2})
    schmidtnumbers = eigvals(dmat)
    vne = 0
    for i in length(schmidtnumbers)
        vne = vne - schmidtnumbers[i]*log(schmidtnumbers[i])
    end
    return vne
end

#Define projection, ladder, identity#
X = [0. 1 0 ; 1 0 1 ; 0 1 0]*(1/sqrt(2))*(1+0*im)
Y = [0. -im 0 ; im 0 -im ; 0 im 0]*(1/sqrt(2))*(1+0*im)
Z = [1. 0 0 ; 0 0 0; 0 0 -1]*(1+0*im)
plus = [1. 0 0 ; 0 0 0; 0 0 0]*(1+0*im)
minus = [0. 0 0 ; 0 1 0; 0 0 0]*(1+0*im)
projplus = [1. 0 0 ; 0 0 0 ; 0 0 0]*(1+0*im)
projzero = [0. 0 0 ; 0 1 0 ; 0 0 0]*(1+0*im)
projminus = [0. 0 0 ; 0 0 0 ; 0 0 1]*(1+0*im)
id = [1. 0 0 ; 0 1 0 ; 0 0 1]*(1+0*im)

#system length#
L=5
#which excited state-- 1-3^L are permitted
n=1


H = zeros(3^L,3^L)
for i in 1:L-1
    H = H + MakeOp(MakeOpList(X,X,i,i+1,L))+MakeOp(MakeOpList(Y,Y,i,i+1,L))+MakeOp(MakeOpList(Z,Z,i,i+1,L))+(1/3)*(MakeOp(MakeOpList(X,X,i,i+1,L))+MakeOp(MakeOpList(Y,Y,i,i+1,L))+MakeOp(MakeOpList(Z,Z,i,i+1,L)))*(MakeOp(MakeOpList(X,X,i,i+1,L))+MakeOp(MakeOpList(Y,Y,i,i+1,L))+MakeOp(MakeOpList(Z,Z,i,i+1,L)))
end

states = eigvecs(H)
dmat = makerho(states,n)
YAX = zeros(L-1)
for i in 1:L-1
    YAX[i] = SvN(partialtrace(dmat,i))
end
XAX = 2:L

plot(XAX,YAX,title = string("SvN of AKLT Chain of Length ",L, " N =",n), label =  nothing, xlabel = "Bipartition Placement", ylabel = "Entanglement Entropy")


#savefig(string("Density of States with L= ",L,", W=", w, ", and t=",t))


#since we assume occupation of 1 per site, we have filled the first L eigenstates#
