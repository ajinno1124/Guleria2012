include("SkyrmeParams.jl")
include("MyLib.jl")
using Parameters
using LinearAlgebra
using SparseArrays
using Arpack
using .NuclParameters
using .LambdaParameters
using .MyLib

@consts begin
    # Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
    mpMeV=938.27208816
    mnMeV=939.5654205
    mAveMeV=(mpMeV+mnMeV)/2
    mΛMeV=1115.683 #\pm 0.006 MeV
    ħc=197.3269804
    e2MeVfm=1.4400 

    Nmesh=500
    rmax=25
    nmax=5
    lmax=7

    offset=100 #Sparse Matrixを計算する際のoffset
end

mutable struct QuantumNumber
    j::Float64
    l::Int64
    B::Int64 #1::p, 2::n, 3::Λ
end

mutable struct SingleParticleState
    QN::QuantumNumber
    E::Float64
    ψ::Vector{Float64}
end

mutable struct AtomNum
    Z::Int64
    N::Int64
    Λ::Int64
end

function getrmesh()
    #h=rmax/(Nmesh-0.5)
    #rmesh=range(0.5*h, (Nmesh-0.5)*h, length=Nmesh)
    rmesh=range(1e-3,rmax,length=Nmesh)
    return rmesh
end

##########################################################
# make Hamiltonian matrix
function Calc_ρ(occ,States::Vector{SingleParticleState},rmesh)
    ρ=zeros(Float64,Nmesh)
    for i in 1:length(occ)
        j=States[i].QN.j
        @. ρ[:]+=occ[i]*(2*j+1)/(4*π)*(States[i].ψ[:]/rmesh[:])
    end
    return ρ
end

function Calc_dρ(ρ,rmesh)
    h=rmesh[2]-rmesh[1]
    dρ=MyLib.diff1st(h,ρ)
    return dρ
end

function Calc_Lapρ(occ,States::Vector{SingleParticleState},rmesh)
    Lapρ=zeros(Float64,Nmesh)
    ρ=Calc_ρ(occ,States,rmesh)
    h=rmesh[2]-rmesh[1]
    dρ=MyLib.diff1st(h,ρ)
    ddρ=MyLib.diff2nd(h,ρ)
    @. Lapρ=2/rmesh[:]*dρ + ddρ
    return Lapρ
end

function Calc_τ(occ,States::Vector{SingleParticleState},rmesh)
    τ=zeros(Float64,Nmesh)
    Rr=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    for i=eachindex(occ)
        @. Rr[:]=States[i].ψ[:]/rmesh[:]
        dRrdr=MyLib.diff1st(h,Rr)
        j=States[i].QN.j
        l=States[i].QN.l
        @. τ[:]+=occ[i]*(2*j+1)/(4*π)*(dRrdr[:]^2 + l*(l+1)/(rmesh[:]^2)*Rr[:]^2)
    end
    return τ
end

function Calc_J(occ,States::Vector{SingleParticleState},rmesh)
    J=zeros(Float64,Nmesh)
    for i=eachindex(occ)
        j=States[i].QN.j
        l=States[i].QN.length
        @. J[:]+=occ[i]/rmesh[i]*(2*j+1)/(4*π)*(j*(j+1)-l*(l+1)-0.75)*(States[i].ψ[:]/rmesh[i])
    end
    return J
end

function Calc_divJ(J,rmesh)
    divJ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    dJdr=MyLib.diff1st(h,J)
    @. divJ+=dJdr[:]+2/rmesh[:]*J[:]
    return divJ
end
function Calc_h2mΛ(aΛ,ρN)
    return @. ħc^2/(2*mΛMeV)+aΛ[2]*ρN
end

function Calc_h2mN(B,aN,aΛ,ρN,ρq,ρΛ)
    m=0.0
    if B==1 #proton
        m+=mnMeV
    elseif B==2 #neutron
        m+=mpMeV
    end
    return @. ħc^2/(2*m)+aN[5]*ρN+aN[6]*ρq+aΛ[2]*ρΛ
end

function Calc_VΛΛ(aΛ,γ,ρN,LapρN,τN)
    return @. aΛ[1]*ρN+aΛ[2]*τN-aΛ[3]*LapρN+aΛ[4]*ρN^(γ+1)
end

function Calc_VΛN(aΛ,γ,ρN,ρΛ,LapρΛ,τΛ,)
    return @. aΛ[1]*ρΛ+aΛ[2]*τΛ-aΛ[3]*LapρΛ+(γ+1)*aΛ[4]*(ρN^γ)*ρΛ
end

function Calc_VNq(aN,σ,W0,ρN,ρq,τN,τq,LapρN,Lapρq,divJN,divJq)
    ans=zeros(Float64,Nmesh)
    @. ans+=2*aN[1]*ρN + 2*aN[2]*ρq + (σ+2)*aN[3]*ρN^(σ+1)
    @. ans+=aN[4]*(σ*ρN^(σ-1)*(ρq^2+(ρN-ρq)^2)+2*(ρN^σ)*ρq)
    @. ans+=aN[5]*τN+aN[6]*τq-2*aN[7]*LapρN
    @. ans+=-2*aN[8]*Lapρq-0.5*W0*(divJN+divJq)
    return ans
end

function Calc_Vcoul(ρp,rmesh,Z)
    Vcoul=zeros(Float64,Nmesh)
    Vcoul+=MyLib.SolvePoissonEq(ρp,rmesh,Z)
    @. Vcoul[:]+=-(3*ρp[:]/π)^(1/3)
    Voul*=e2MeVfm

    return Vcoul
    
end

function Calc_Wq(aN,dρN,dρq,JN,Jq)
    return @. 0.5*(dρN+dρq)+2*aN[9]*JN+2*aN[10]*Jq
end

#Ay"+By'+Cy=ϵy
#=
function CoefA(h2mVal)
    return -h2mVal
end

function CoefB(dh2mVal)
    return dh2mVal
end

function CoefC(h2mVal,dh2mVal,Vq,Wq,QN::QuantumNumber)
    ans=zeros(Float64,Nmesh)
    rmesh=getrmesh()
    @. ans+=h2mVal*QN.l*(QN.l+1)/rmesh
    @. ans+=Vq
    @. ans+=Wq/rmesh*(QN.j*(QN.j+1)-QN.l*(QN.l+1)-0.75)
    return ans
end
=#

function calc_cij2(i,j,h)
    cij = 0.0
    if i==1
        cij += ifelse(j==1,-3.0,0.0)
        cij += ifelse(j==2,4.0,0.0)
        cij += ifelse(j==3,-1.0,0.0)
    elseif i==Nmesh
        cij += ifelse(j==Nmesh,3.0,0.0)
        cij += ifelse(j==Nmesh-1,-4.0,0.0)
        cij += ifelse(j==Nmesh-2,1.0,0.0)
    else
        cij += ifelse(j==i+1,1.0,0.0)
        cij += ifelse(j==i-1,-1.0,0.0)
    end
    cij = cij/(2*h)
    return cij
end

function calc_dij2(i,j,h)
    dij = 0.0
    if i==1
        dij += ifelse(j==1,2.0,0.0)
        dij += ifelse(j==2,-5.0,0.0)
        dij += ifelse(j==3,4.0,0.0)
        dij += ifelse(j==4,-1.0,0.0)
    elseif i==Nmesh
        dij += ifelse(j==Nmesh,2.0,0.0)
        dij += ifelse(j==Nmesh-1,-5.0,0.0)
        dij += ifelse(j==Nmesh-2,4.0,0.0)
        dij += ifelse(j==Nmesh-3,-1.0,0.0)
    else
        dij += ifelse(j==i+1,1.0,0.0)
        dij += ifelse(j==i,-2.0,0.0)
        dij += ifelse(j==i-1,1.0,0.0)
    end
    dij = dij/(h^2)
    return dij
end

#functionに入れる段階でBは決めておく
function MakeHmat(A,B,C,rmesh)
    Hmat=spzeros(Nmesh,Nmesh)
    h=rmesh[2]-rmesh[1]

    for row in 1:Nmesh
        Hmat[row,row]+=C[row]
        for col in max(row-3,1):min(row+3,Nmesh)
            cij=calc_cij2(row,col,h)
            dij=calc_dij2(row,col,h)
            if cij!=0.0 || dij!=0.0
                Hmat[row,col]+=A[row]*dij+B[row]*cij
            end
        end

    end

    return Hmat
end

##########################################################
# HF iteration

function NormFact(rmesh,ψ)
    ans=sqrt(MyLib.IntTrap(rmesh,@. ψ[:]^2))
    return 1/ans
end

function CalcStates(QN::QuantumNumber,h2mB,dh2mB,VB,WB,rmesh)
    States=SingleParticleState[]

    j=QN.j
    l=QN.l
    A=zeros(Float64,Nmesh)
    B=zeros(Float64,Nmesh)
    C=zeros(Float64,Nmesh)
    @. A[:] += -h2mB[:]
    @. B[:] += -dh2mB[:]
    @. C[:] += h2mB[:]*l*(l+1)/(rmesh[:]^2) + VB[:] + dh2mB[:]/rmesh[:] + WB[:]/rmesh[:]*(j*(j+1)-l*(l+1)-0.75)


    Hmat=MakeHmat(A,B,C,rmesh)
    for i in 1:Nmesh
        Hmat[i,i]+=offset
    end

    E,ψ=eigs(Hmat,nev=nmax,which=:SR,maxiter=3000)
    
    E=real(E.-offset)
    for i=eachindex(E)
        Norm=NormFact(rmesh,real(ψ[:,i]))
        @. ψ[:,i]=real(ψ[:,i])/Norm
        if i>1 && E[i-1]!=E[i]
            push!(States,SingleParticleState(QN,real(E[i]),ψ[:,i]))
        end
    end

    return States
end

function CalcAllStates(h2m,dh2m,V,W,rmesh)
    AllStates=[SingleParticleState[],SingleParticleState[],SingleParticleState[]]
    for b in 1:3
        States=SingleParticleState[]
        for l in 0:lmax
            for j in max(l-0.5,0.5):l+0.5
                QN=QuantumNumber(j,l,b)
                States=vcat(States,CalcStates(QN,h2m[b,:],dh2m[b,:],V[b,:],W[b,:],rmesh))
            end
        end
        sort!(States, by=x->x.E)
        AllStates[b]=vcat(AllStates[b],States)
    end

    return AllStates
end

function Calc_occ(AN::AtomNum,AllStates;BCS=false)
    Remain=AtomNum
    Allocc=[Float64[],Float64[],Float64[]]
    if BCS==false
        for b in 1:3
            for i in 1:length(AllStates[b])
                j=AllStates[b][i].QN.j
                if Remain[b]>=(2*j+1)
                    push!(occ[b],1)
                    Remain[b]-=2*j+1
                else
                    push!(occ[b],Remain[b]/(2*j+1))
                    Remain[b]=0
                end
            end
        end
    end

    return Allocc
end

function InitPot(AN::AtomNum,rmesh)
    N=AN.N
    Z=AN.Z
    A=N+Z
    Vcoef=zeros(Float64,3)
    Vcoef[1]=-52-33*((N-Z)/A) #proton
    Vcoef[2]=-52+33*((N-Z)/A) #neutron
    Vcoef[3]=-52*(2/3) #Lambda

    r0=1.27
    a=0.67
    R=r0*A^(1/3)
    h2m=zeros(Float64,(3,Nmesh))
    dh2m=zeros(Float64,(3,Nmesh))
    V=zeros(Float64,(3,Nmesh))
    W=zeros(Float64,(3,Nmesh))

    for b in 1:3
        if b==1 #proton
            @. h2m[b,:]+=ħc^2/(2*mpMeV)
            @. V[b,:]+=Vcoef[b]/(1+exp((rmesh[:]-R)/a))
            @. W[b,:]+=(0.22*Vcoef[b,:]*r0^2/a)*(exp((rmesh[:]-R)/a)/(1+exp((rmesh[:]-R)/a))^2)
        elseif b==2
            @. h2m[b,:]+=ħc^2/(2*mnMeV)
            @. V[b,:]+=Vcoef[b]/(1+exp((rmesh[:]-R)/a))
            @. W[b,:]+=(0.22*Vcoef[b,:]*r0^2/a)*(exp((rmesh[:]-R)/a)/(1+exp((rmesh[:]-R)/a))^2)
        elseif b==3
            @. h2m[b,:]+=ħc^2/(2*mΛMeV)
            @. V[b,:]+=Vcoef[b]/(1+exp((rmesh[:]-R)/a))
            @. W[b,:]+=0
        end
    end

    return h2m,dh2m,V,W
end

# create initial state using WoodsSaxon Potential
function InitialCondition(AN::AtomNum)
    rmesh=getrmesh()
    h2m,dh2m,V,W=InitPot(AN,rmesh)

    InitState=CalcAllStates(h2m,dh2m,V,W,rmesh)
    #Initocc=Calc_occ(AN,InitState)

    return InitState
end

function CheckConvergence(OldStates::Vector{SingleParticleState},NewStates::Vector{SingleParticleState};rtol=1e-4)
    
end

function Calc_Density(Allocc,AllStates)
    ρ3=zeros(Float64,(3,Nmesh))
    dρ3=zeros(Float64,(3,Nmesh))
    Lapρ3=zeros(Float64,(3,Nmesh))
    τ3=zeros(Float64,(3,Nmesh))
    J3=zeros(Float64,(3,Nmesh))
    divJ3=zeros(Float64,(3,Nmesh))
    rmesh=getrmesh()

    for b in 1:3
        ρ3[b,:]=Calc_ρ(Allocc[b],AllStates[b],rmesh)
        dρ3[b,:]=Calc_dρ(ρ3[b,:],rmesh)
        Lapρ3[b,:]=Calc_Lapρ(Allocc[b],AllStates[b],rmesh)
        τ3[b,:]=Calc_τ(Allocc[b],AllStates[b],rmesh)
        J3[b,:]=Calc_J(Allocc[b],AllStates[b],rmesh)
        divJ3[b,:]=Calc_divJ(J,rmesh)
    end

    return ρ3,dρ3,Lapρ3,τ3,J3,divJ3
end

function Calc_Coef(occ,AllStates,aN,aΛ,pN,pΛ,Z)
    rmesh=getrmesh()

    ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(occ,AllStates)
    ρN=ρ3[1]+ρ3[2]
    dρN=dρ3[1]+dρ3[2]
    LapρN=Lapρ3[1]+Lapρ3[2]
    JN=J3[1]+J[2]
    divJN=divJ3[1]+divJ3[2]
    h=rmesh[2]-rmesh[1]

    h2m=zeros(Float64,(3,Nmesh))
    dh2m=zeros(Float64,(3,Nmesh))
    V=zeros(Float64,(3,Nmesh))
    W=zeros(Float64,(3,Nmesh))

    VΛΛ=Calc_VΛΛ(aΛ, pΛ.γ, ρN,LapρN,τN)
    VΛN=Calc_VΛN(aΛ, pΛ.γ, ρN, ρ3[3],Lapρ3[3],τ3[3])
    VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1], τN, τ3[1],LapρN,Lapρ3[1],divJN,divJ3[1])
    VNn=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[2], τN, τ3[2],LapρN,Lapρ3[2],divJN,divJ3[2])
    Vcoul=Calc_Vcoul(ρp,rmesh,Z)

    for b in 1:3
        if b==1 #proton
            h2m[b,:]+=Calc_h2mN(b,aN,aΛ,ρN,ρ3[b],ρ3[3])
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛN+VNp+Vcoul
            W[b,:]+=Calc_Wq(aN,dρN,dρ3[b],JN,J3[b])
        elseif b==2 #neutron
            h2m[b,:]+=Calc_h2mN(b,aN,aΛ,ρN,ρ3[b],ρ3[3])
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛN+VNn
            W[b,:]+=Calc_Wq(aN,dρN,dρ3[b],JN,J3[b])
        elseif b==3 #Lambda
            h2m[b,:]+=Calc_h2mΛ(aΛ,ρN)
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛΛ
            W[b,:]+=0
        end
    end

    return h2m,dh2m,V,W
end

function HF_iter(AN::AtomNum;MaxIter=10,NParamType="SLy4",ΛParamType="HPΛ1")
    OldStates=InitialCondition(AN)
    Oldocc=Calc_occ(AN,OldStates)
    rmesh=getrmesh()

    aN=NuclParameters.getaN(NParamType)
    aΛ=LambdaParameters.getaΛ(ΛParamType)
    pN=NuclParameters.getParams(NParamType)
    pΛ=LambdaParameters.getParams(ΛParamType)

    #calc several params
    for i in 1:MaxIter
        h2m,dh2m,V,W=Calc_Coef(Oldocc,AllStates,aN,aΛ,pN,pΛ,Z)
        NewStates=CalcAllStates(h2m,dh2m,V,W,rmesh)
        Newocc=Calc_occ(AN,NewStates)
        
        if CheckConvergence(OldStates,NewStates)==true
            break
        end

        OldStates=NewStates
        Oldocc=Newocc
    end

    return NewStates
end



