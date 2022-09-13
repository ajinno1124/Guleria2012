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

    Nmesh=5000
    Nmatch=2000
    rmax=30
    nmax=4
    lmax=6

    offset=1000 #Sparse Matrixを計算する際のoffset
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
    rmesh=range(1e-2,rmax,length=Nmesh)
    return rmesh
end

##########################################################
# Defenition of Density, Potential
function Calc_ρ(occ::Vector{Float64},States::Vector{SingleParticleState},rmesh)
    ρ=zeros(Float64,Nmesh)
    for i in 1:length(occ)
        j=States[i].QN.j
        @. ρ[:]+=occ[i]*(2*j+1)/(4*π)*(States[i].ψ[:]/rmesh[:])^2
    end
    return ρ
end

function Calc_dρ(ρ,rmesh)
    h=rmesh[2]-rmesh[1]
    dρ=MyLib.diff1st(h,ρ)
    return dρ
end

function Calc_Lapρ(ρ::Vector{Float64},rmesh)
    Lapρ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    dρ=MyLib.diff1st(h,ρ)
    ddρ=MyLib.diff2nd(h,ρ)
    @. Lapρ[:]=2/rmesh[:]*dρ[:] + ddρ[:]
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
        l=States[i].QN.l
        @. J[:]+=occ[i]/rmesh[i]*(2*j+1)/(4*π)*(j*(j+1)-l*(l+1)-0.75)*(States[i].ψ[:]/rmesh[:])^2
    end
    return J
end

function Calc_divJ(J::Vector{Float64},rmesh)
    divJ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    dJdr=MyLib.diff1st(h,J)
    @. divJ+=dJdr[:]+2/rmesh[:]*J[:]
    return divJ
end
function Calc_h2mΛ(aΛ,ρN::Vector{Float64})
    return @. ħc^2/(2*mΛMeV)+aΛ[2]*ρN
end

function Calc_h2mN(B,aN,aΛ,ρN::Vector{Float64},ρq::Vector{Float64},ρΛ::Vector{Float64})
    m=0.0
    if B==1 #proton
        m+=mnMeV
    elseif B==2 #neutron
        m+=mpMeV
    end
    return @. ħc^2/(2*m)+aN[5]*ρN[:]+aN[6]*ρq[:]+aΛ[2]*ρΛ[:]
end

function Calc_VΛΛ(aΛ,γ,ρN::Vector{Float64},LapρN::Vector{Float64},τN::Vector{Float64})
    return @. aΛ[1]*ρN+aΛ[2]*τN-aΛ[3]*LapρN+aΛ[4]*ρN^(γ+1)
end

function Calc_VΛN(aΛ,γ,ρN::Vector{Float64},ρΛ::Vector{Float64},LapρΛ::Vector{Float64},τΛ::Vector{Float64})
    return @. aΛ[1]*ρΛ+aΛ[2]*τΛ-aΛ[3]*LapρΛ+(γ+1)*aΛ[4]*(ρN^γ)*ρΛ
end

function Calc_VNq(aN,σ,W0,ρN,ρq,τN,τq,LapρN,Lapρq,divJN,divJq)
    ans=zeros(Float64,Nmesh)
    @. ans[:]+=2*aN[1]*ρN[:] + 2*aN[2]*ρq[:] + (σ+2)*aN[3]*ρN[:]^(σ+1)
    @. ans[:]+=aN[4]*(σ*ρN[:]^(σ-1)*(ρq[:]^2+(ρN[:]-ρq[:])^2)+2*(ρN[:]^σ)*ρq[:])
    @. ans[:]+=aN[5]*τN[:]+aN[6]*τq[:]-2*aN[7]*LapρN[:]
    @. ans[:]+=-2*aN[8]*Lapρq[:]-0.5*W0*(divJN[:]+divJq[:])
    return ans
end

function Calc_Vcoul(ρp::Vector{Float64},rmesh,Z)
    Vcoul=zeros(Float64,Nmesh)
    Vcoul+=MyLib.SolvePoissonEq(ρp,rmesh,Z)
    @. Vcoul[:]=Vcoul[:]/rmesh[:]
    @. Vcoul[:]+=-(3*ρp[:]/π)^(1/3)
    Vcoul*=e2MeVfm

    return Vcoul
    
end

function Calc_Wq(aN,W0,dρN::Vector{Float64},dρq::Vector{Float64},JN::Vector{Float64},Jq::Vector{Float64})
    return @. 0.5*W0*(dρN[:]+dρq[:])+2*aN[9]*JN[:]+2*aN[10]*Jq[:]
end


##########################################################
# Euler Method
function Euler2(y,A,B,h)
	ynext=(1-h*A)*y+h*B
	return ynext
end

#v'=u
function DiffEqOfv(v,u,h)
	A=
	B=u
	return Euler2(v,A,B,h)
end

#
function DiffEqOfu(u,v,E,V,r,l,h)
    A=
	B=-(2*mAveMeV*(E-V)/ħc^2-l*(l+1)/r^2)*v
	return Euler2(u,A,B,h)
end

function InitCondEuler(l,E,rmesh)
    h=rmesh[2]-rmesh[1]
    vin1=rmesh[1]^(l+1)
    uin1=(l+1)*rmesh[1]^l
    vout2=exp(-(-2*mAveMeV/ħc^2*E)^0.5*rmesh[Nmesh])
    uout2=-(-2*mAveMeV/ħc^2*E)^0.5*exp(-(-2*mAveMeV/ħc^2*E)^0.5*rmesh[Nmesh])
    return vin1,uin1,vout2,uout2 #あとでnormalizeする。
end

function ShootingEuler(l,E,V::Vector{Float64},rmesh)
    vin=zeros(Float64,2)
    vout=zeros(Float64,2)
    uin=zeros(Float64,2)
    uout=zeros(Float64,2)

    rmesh=getrmesh()
    h=rmesh[2]-rmesh[1]
    vin[1],uin[1],vout[2],uout[2]=InitCondEuler(l,E,rmesh)

    for i in 1:Nmatch-1
        vin[2]=DiffEqOfv(vin[1],uin[1],h)
        uin[2]=DiffEqOfu(uin[1],vin[1],E,V[i],rmesh[i],l,h)
        vin[1]=vin[2]
        uin[1]=uin[2]
    end

    for i in Nmesh:-1:Nmatch-1
        vout[1]=DiffEqOfv(vout[2],uout[2],-h)
        uout[1]=DiffEqOfu(uout[2],vout[2],E,V[i],rmesh[i],l,-h)
        vout[2]=vout[1]
        uout[2]=uout[1]
    end

    return vin[1]*uout[1]-vout[1]*uin[1]
    
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

    
    for i in 1:nmax
        E,ψ=
        Norm=NormFact(rmesh,real(ψ[:,i]))
        @. ψ[:,i]=sign(ψ[2,i]-ψ[1,i])*Norm*real(ψ[:,i])
        push!(States,SingleParticleState(QN,real(E[i]),ψ[:,i]))
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
    Allocc=[Float64[],Float64[],Float64[]]
    Remain=[AN.Z,AN.N,AN.Λ]
    if BCS==false
        for b in 1:3
            for i=eachindex(AllStates[b])
                j=AllStates[b][i].QN.j
                if Remain[b]>=(2*j+1)
                    push!(Allocc[b],1.0)
                    Remain[b]-=2*j+1
                else
                    push!(Allocc[b],Remain[b]/(2*j+1))
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
            @. W[b,:]+=(0.22*Vcoef[b]*r0^2/a)*(exp((rmesh[:]-R)/a)/(1+exp((rmesh[:]-R)/a))^2)
        elseif b==2
            @. h2m[b,:]+=ħc^2/(2*mnMeV)
            @. V[b,:]+=Vcoef[b]/(1+exp((rmesh[:]-R)/a))
            @. W[b,:]+=(0.22*Vcoef[b]*r0^2/a)*(exp((rmesh[:]-R)/a)/(1+exp((rmesh[:]-R)/a))^2)
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

function CheckConvergence(Oldocc,OldStates,Newocc,NewStates,rmesh;rtol=1e-2)
    check=true
    diffρ=zeros(Float64,3)
    for b in 1:3
        ρold=Calc_ρ(Oldocc[b],OldStates[b],rmesh)
        ρnew=Calc_ρ(Newocc[b],NewStates[b],rmesh)
        #diffρ[b]=MyLib.IntTrap(rmesh,@. (4*π*rmesh[:]^2*(ρold[:]-ρnew[:])))
        diffρ[b]=MyLib.IntTrap(rmesh,(@. ((ρold[:]-ρnew[:])/(ρold[:]+ρnew[:]))^2))
        if abs(diffρ[b])>rtol
            check=false
        end
    end

    println("diffρ=($(diffρ[1]),$(diffρ[2]),$(diffρ[3]))")

    return check

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
        Lapρ3[b,:]=Calc_Lapρ(ρ3[b,:],rmesh)
        τ3[b,:]=Calc_τ(Allocc[b],AllStates[b],rmesh)
        J3[b,:]=Calc_J(Allocc[b],AllStates[b],rmesh)
        divJ3[b,:]=Calc_divJ(J3[b,:],rmesh)
    end

    return ρ3,dρ3,Lapρ3,τ3,J3,divJ3
end

function Calc_Density2(Allocc,AllStates)
    ρ3=zeros(Float64,(3,Nmesh))
    τ3=zeros(Float64,(3,Nmesh))
    J3=zeros(Float64,(3,Nmesh))
    rmesh=getrmesh()

    for b in 1:3
        ρ3[b,:]=Calc_ρ(Allocc[b],AllStates[b],rmesh)
        τ3[b,:]=Calc_τ(Allocc[b],AllStates[b],rmesh)
        J3[b,:]=Calc_J(Allocc[b],AllStates[b],rmesh)
    end

    return ρ3,τ3,J3
end

function Calc_Coef(occ,AllStates,aN,aΛ,pN,pΛ,Z)
    rmesh=getrmesh()

    ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(occ,AllStates)
    ρN=ρ3[1,:]+ρ3[2,:]
    dρN=dρ3[1,:]+dρ3[2,:]
    LapρN=Lapρ3[1,:]+Lapρ3[2,:]
    τN=τ3[1,:]+τ3[2,:]
    JN=J3[1,:]+J3[2,:]
    divJN=divJ3[1,:]+divJ3[2,:]
    h=rmesh[2]-rmesh[1]

    h2m=zeros(Float64,(3,Nmesh))
    dh2m=zeros(Float64,(3,Nmesh))
    V=zeros(Float64,(3,Nmesh))
    W=zeros(Float64,(3,Nmesh))

    VΛΛ=Calc_VΛΛ(aΛ, pΛ.γ, ρN,LapρN,τN)
    VΛN=Calc_VΛN(aΛ, pΛ.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:])
    VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1,:], τN, τ3[1,:],LapρN,Lapρ3[1,:],divJN,divJ3[1,:])
    VNn=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[2,:], τN, τ3[2,:],LapρN,Lapρ3[2,:],divJN,divJ3[2,:])
    Vcoul=Calc_Vcoul(ρ3[1,:],rmesh,Z)

    for b in 1:3
        if b==1 #proton
            h2m[b,:]+=Calc_h2mN(b,aN,aΛ,ρN,ρ3[b,:],ρ3[3,:])
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛN+VNp+Vcoul
            W[b,:]+=Calc_Wq(aN,pN.W0,dρN,dρ3[b,:],JN,J3[b,:])
        elseif b==2 #neutron
            h2m[b,:]+=Calc_h2mN(b,aN,aΛ,ρN,ρ3[b,:],ρ3[3,:])
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛN+VNn
            W[b,:]+=Calc_Wq(aN,pN.W0,dρN,dρ3[b,:],JN,J3[b,:])
        elseif b==3 #Lambda
            h2m[b,:]+=Calc_h2mΛ(aΛ,ρN)
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛΛ
            @. W[b,:]+=0
        end
    end

    return h2m,dh2m,V,W
end

function Calc_Coef(ρ3,τ3,J3,aN,aΛ,pN,pΛ,Z)
    dρ3=zeros(Float64,(3,Nmesh))
    Lapρ3=zeros(Float64,(3,Nmesh))
    divJ3=zeros(Float64,(3,Nmesh))
    rmesh=getrmesh()

    for b in 1:3
        dρ3[b,:]=Calc_dρ(ρ3[b,:],rmesh)
        Lapρ3[b,:]=Calc_Lapρ(ρ3[b,:],rmesh)
        divJ3[b,:]=Calc_divJ(J3[b,:],rmesh)
    end

    ρN=ρ3[1,:]+ρ3[2,:]
    dρN=dρ3[1,:]+dρ3[2,:]
    LapρN=Lapρ3[1,:]+Lapρ3[2,:]
    τN=τ3[1,:]+τ3[2,:]
    JN=J3[1,:]+J3[2,:]
    divJN=divJ3[1,:]+divJ3[2,:]
    h=rmesh[2]-rmesh[1]

    h2m=zeros(Float64,(3,Nmesh))
    dh2m=zeros(Float64,(3,Nmesh))
    V=zeros(Float64,(3,Nmesh))
    W=zeros(Float64,(3,Nmesh))

    VΛΛ=Calc_VΛΛ(aΛ, pΛ.γ, ρN,LapρN,τN)
    VΛN=Calc_VΛN(aΛ, pΛ.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:])
    VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1,:], τN, τ3[1,:],LapρN,Lapρ3[1,:],divJN,divJ3[1,:])
    VNn=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[2,:], τN, τ3[2,:],LapρN,Lapρ3[2,:],divJN,divJ3[2,:])
    Vcoul=Calc_Vcoul(ρ3[1,:],rmesh,Z)

    for b in 1:3
        if b==1 #proton
            h2m[b,:]+=Calc_h2mN(b,aN,aΛ,ρN,ρ3[b,:],ρ3[3,:])
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛN+VNp+Vcoul
            W[b,:]+=Calc_Wq(aN,pN.W0,dρN,dρ3[b,:],JN,J3[b,:])
        elseif b==2 #neutron
            h2m[b,:]+=Calc_h2mN(b,aN,aΛ,ρN,ρ3[b,:],ρ3[3,:])
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛN+VNn
            W[b,:]+=Calc_Wq(aN,pN.W0,dρN,dρ3[b,:],JN,J3[b,:])
        elseif b==3 #Lambda
            h2m[b,:]+=Calc_h2mΛ(aΛ,ρN)
            dh2m[b,:]+=MyLib.diff1st(h,h2m[b,:])
            V[b,:]+=VΛΛ
            @. W[b,:]+=0
        end
    end

    return h2m,dh2m,V,W
end

function HF_iter(AN::AtomNum;MaxIter=60,NParamType="SLy4",ΛParamType="HPΛ1")
    OldStates=InitialCondition(AN)
    Oldocc=Calc_occ(AN,OldStates)
    rmesh=getrmesh()
    Oldρ3,Oldτ3,OldJ3=Calc_Density2(Oldocc,OldStates)

    aN=NuclParameters.getaN(NParamType)
    aΛ=LambdaParameters.getaΛ(ΛParamType)
    pN=NuclParameters.getParams(NParamType)
    pΛ=LambdaParameters.getParams(ΛParamType)

    #for debug
    ρptest=zeros(Float64,(MaxIter,Nmesh))

    #calc several params
    for i in 1:MaxIter
        #for debug
        ρptest[i,:]=Calc_ρ(Oldocc[1],OldStates[1],rmesh)

        h2m,dh2m,V,W=Calc_Coef(Oldρ3,Oldτ3,OldJ3,aN,aΛ,pN,pΛ,AN.Z)

        NewStates=CalcAllStates(h2m,dh2m,V,W,rmesh)
        Newocc=Calc_occ(AN,NewStates)
        Newρ3,Newτ3,NewJ3=Calc_Density2(Newocc,NewStates)
        
        if CheckConvergence(Oldocc,OldStates,Newocc,NewStates,rmesh)==true
            return Newocc,NewStates
            #break
        end
        println(i)
        OldStates=NewStates
        Oldocc=Newocc
        α=0.1
        Oldρ3=Oldρ3*(1-α)+Newρ3*α
        Oldτ3=Oldτ3*(1-α)+Newτ3*α
        OldJ3=OldJ3*(1-α)+NewJ3*α
    end

    plot(xlabel="r",ylabel="ρn")
    for i in 1:MaxIter
        plot!(rmesh,ρptest[i,:],label="$i")
    end
    plot!()

end



