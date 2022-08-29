include("SkyrmeParams.jl")
using Parameters
using LinearAlgebra
using Arpack

@consts begin
    # Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
    mpMeV=938.27208816
    mnMeV=939.5654205
    mAveMeV=(mpMeV+mnMeV)/2
    mΛMeV=1115.683 #\pm 0.006 MeV
    ħc=197.3269804

    Nmesh=1000
    rmax=15
    nmax=10
    lmax=7

    offset=100 #Sparse Matrixを計算する際のoffset
end

mutable struct QuantumNumber
    j::Float64
    l::Int64
    τ::Int64 #-1::p, 0::Lambda, 1::Ptoron
end

mutable struct SingleParticleState
    QN::QuantumNumber
    E::Float64
    ψ::Vector{Float64}
end

function getrmesh()
    h=rmax/(Nmesh-0.5)
    rmesh=0.5*h:h:(Nmesh-0.5)*h
    return rmesh
end

#Ay"+By'+Cy=ϵy
function CoefA()
end

function CoefB()
end

function CoefC()
end

function calc_cij2(i,j,h,P)
    cij = 0.0
    if i==1
        cij += ifelse(j==2,1.0,0.0)        
        cij += ifelse(j==1,-P,0.0)
    elseif i==Nmesh
        #cij += ifelse(j==Nmesh-1,-1.0,0.0)
        cij += ifelse(j==Nmesh-1,-2.0,0.0)
        cij += ifelse(j==Nmesh,2.0,0.0)
    else
        cij += ifelse(j==i+1,1.0,0.0)
        cij += ifelse(j==i-1,-1.0,0.0)
    end
    cij = cij/(2*h)
    return cij
end

function calc_dij4(i,j,h)
    dij = 0.0
    if i==1
        dij += ifelse(j==2,1.0,0.0)
        dij += ifelse(j==1,P-2.0,0.0)
    elseif i==Nmesh
        #dij += ifelse(j==Nmesh,-2.0,0.0)        
        #dij += ifelse(j==Nmesh-1,1.0,0.0)

        dij += ifelse(j==Nmesh,2.0,0.0)
        dij += ifelse(j==Nmesh-1,-5.0,0.0)
        dij += ifelse(j==Nmesh-2,4.0,0.0)
        dij += ifelse(j==Nmesh-3,-1.0,0.0)

        #dij += ifelse(j==Nmesh,1.0,0.0)
        #dij += ifelse(j==Nmesh-1,-2.0,0.0)
        #dij += ifelse(j==Nmesh-2,1.0,0.0)
    else
        dij += ifelse(j==i+1,1.0,0.0)
        dij += ifelse(j==i,-2.0,0.0)
        dij += ifelse(j==i-1,1.0,0.0)
    end
    dij = dij/(h^2)
    return dij
end

function MakeHmat(rmesh)
    Hmat=spzeros(Nmesh,Nmesh)
    h=rmesh[2]-rmesh[1]

    for i in 1:Nmesh
        Ai=CoefA()
        Bi=CoefB()
        Ci=CoefC()
        Hmat[i,i]+=Ci
        
        for j in max(i-3,1):min(i+3,Nmesh)
            cij=calc_cij2(i,j,h,l)
            dij=calc_dij2(i,j,h,l)
            if cij!=0 || dij!=0
                Hmat[i,j]+=Ai*dij+Bi*cij
            end
        end

    end

    return Hmat
end

function IntTrap(x,y)
    N=length(x)
    ans=0.0
    ans+=2*sum(y)-y[N]-y[1]
    ans*=(x[N]-x[1])/(2*(N-1))
    return ans
end

function NormFact(rmesh,ψ)
    ans=sqrt(IntTrap(rmesh,ψ))
    return 1/ans
end

function CalcStates(QN::QuantumNumber)
    States=SingleParticleState[]
    rmesh=getrmesh()

    Hmat=MakeHmat(rmesh)
    for i in 1:Nmesh
        Hmat[i,i]+=offset
    end

    E,ψ=eigs(Hmat,nev=nmax)
    E=real(E.-offset)
    for i in 1:nmax
        ψ[:,i]=real(ψ[:,i])./NormFact(rmesh,real(ψ[:,i]))
        push!(States,SingleParticleState(QN,E[i],ψ[:,i]))
    end

    return States
end

function CalcAllStates(τ)
    AllStates=SingleParticleState[]
    for l in 0:lmax
        for j in max(l-0.5,0.5):l+0.5
            QN=QuantumNumber(j,l,τ)
            States=CalcStates(QN)
            AllStates=vcat(AllStates,States)
        end
    end

    sort!(AllStates, by=x->x.E)

    return AllStates
end

function occ()
end

function InitialCondition()
end

function Calc_ρ()
end

function Calc_τ()
end

function Calc_J()
end

function HF_iter()
    InitialCondition()
end



