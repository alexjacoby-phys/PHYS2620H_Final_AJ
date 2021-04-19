using LinearAlgebra
using Combinatorics
using Plots


#creates list of operators to pass to makeop; dims option = True only if local hilbert space is not uniform-- pending#
function MakeOpList(A::Array{Float64,2},B::Array{Float64,2},x::Int64, y::Int64,L::Int64)
    Id = [1. 0 0 ; 0 1 0; 0 0 1]
    OpList = Vector{Array{Float64,2}}(undef,L)
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





#does tensor product to make operator act on the whole space#
function MakeOp(C::Vector{Array{Float64,2}})
    HTerm = C[1]
    for i in 2:length(C)
        HTerm = kron(HTerm,C[i])
    end
    return HTerm
end


#Define projection, ladder, identity#
X = [0. 1 0 ; 1 0 1 ; 0 1 0]
Y = [0. -im 0 ; im 0 -im ; 0 im 0]
Z = [1. 0 0 ; 0 0 0; 0 0 -1]
plus = [1. 0 0 ; 0 0 0; 0 0 0]
minus = [0. 0 0 ; 0 1 0; 0 0 0]
projplus =
projzero
projminus
id = [1. 0 0 ; 0 1 0 ; 0 0 1]

#system length#
L=10




H = zeros(3^L,3^L)
for i in 1:L-1
    H = H + t*MakeOp(MakeOpList(plus,minus,i,i+1,L))+t*MakeOp(MakeOpList(minus,plus,i,i+1,L))+OnSite[i]*MakeOp(MakeOpList(up,up,i,i,L))
end
H= H + OnSite[L]*MakeOp(MakeOpList(up,up,L,L,L))

Spectrum = eigvals(H)
SpectrumPos = Spectrum - Spectrum[1]*ones(length(Spectrum))






plot(SpectrumPOS,1:length(Spectrum),title = string("Density of States with L= ",L,", W=", w, ", and t=",t), label =  nothing, seriestype = :scatter, xlabel = "Energy", ylabel = "Approximate DOS")


#savefig(string("Density of States with L= ",L,", W=", w, ", and t=",t))


#since we assume occupation of 1 per site, we have filled the first L eigenstates#
