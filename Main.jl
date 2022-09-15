include("SkyrmeParams.jl")
include("MyLib.jl")
using Parameters
using LinearAlgebra
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

    Nmesh=4000
    Nmatch=1200
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
    #Be careful for the SolvePoissionEq
    h=rmax/Nmesh
    rmesh=h:h:h*Nmesh
    return rmesh
end

function getmass(QN::QuantumNumber)
    mass=0.0
    if QN.B==1
        mass+=mpMeV
    elseif QN.B==2
        mass+=mnMeV
    elseif QN.B==3
        mass+=mΛMeV
    end
    return mass
end

#Initail Condition (Woods-Saxon)
function InitPot(AN::AtomNum,rmesh)
    N=AN.N
    Z=AN.Z
    A=N+Z
    Vcoef=zeros(Float64,3)
    Vcoef[1]=-51-33*((N-Z)/A) #proton
    Vcoef[2]=-52+33*((N-Z)/A) #neutron
    Vcoef[3]=-52*(2/3) #Lambda

    r0=1.25
    a=0.65
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

#######################################
# Differential Equation AR''+BR+CR=ER #
# R=u, R'=v                           #
#######################################
function CalcABC(QN,h2mB,dh2mB,VB,WB,rmesh)
    j=QN.j
    l=QN.l
    A=zeros(Float64,Nmesh)
    B=zeros(Float64,Nmesh)
    C=zeros(Float64,Nmesh)
    @. A[:] += -h2mB[:]
    @. B[:] += -dh2mB[:]
    @. C[:] += h2mB[:]*l*(l+1)/(rmesh[:]^2) + VB[:] + dh2mB[:]/rmesh[:] + WB[:]/rmesh[:]*(j*(j+1)-l*(l+1)-0.75)

    return A,B,C
end

# Euler Method
# backward difference
# u'=v
# v' + B/A*v = (E-C)/A*u
function DiffEqBD(u,v,E,A,B,C,h)
    det=1+h*B/A-h^2*(E-C)/A
    unew=((1+h*B/A)*u + h*v)/det
    vnew=(h*(E-C)/A*u + v)/det
    return unew, vnew
end


function BoundCond(QN,E,A1,B1,C1,rmesh)
    l=QN.l
    mass=getmass(QN)
    h=rmesh[2]-rmesh[1]
    u0=0
    v0=1
    uin1,vin1=DiffEqBD(u0,v0,E,A1,B1,C1,h)
    uout2=exp(-(-2*mass/ħc^2*E)^0.5*rmesh[Nmesh])
    vout2=-(-2*mass/ħc^2*E)^0.5*exp(-(-2*mass/ħc^2*E)^0.5*rmesh[Nmesh])
    return uin1,vin1,uout2,vout2
end

function WronskyEuler(E,QN::QuantumNumber,A,B,C,rmesh)
    h=rmesh[2]-rmesh[1]

    uin=zeros(Float64,2)
    vin=zeros(Float64,2)
    uout=zeros(Float64,2)
    vout=zeros(Float64,2)
    uin[1],vin[1],uout[2],vout[2]=BoundCond(QN,E,A[1],B[1],C[1],rmesh)

    for i in 2:Nmatch
        uin[2],vin[2]=DiffEqBD(uin[1],vin[1],E,A[i],B[i],C[i],h)
        uin[1]=uin[2]
        vin[1]=vin[2]
    end

    for i in Nmesh-1:-1:Nmatch
        uout[1],vout[1]=DiffEqBD(uout[2],vout[2],E,A[i],B[i],C[i],-h)
        uout[2]=uout[1]
        vout[2]=vout[1]
    end

    return uin[1]*vout[1]-uout[1]*vin[1]
end

##########################################################
#Calculate States by given A, B, C
function NormFact(rmesh,ψ)
    ans=MyLib.IntTrap(rmesh,@. ψ[:]^2)
    ans+=ψ[1]^2*rmesh[1]/2 #r=0~rmesh[1]までの量を足す
    ans=sqrt(ans)
    return ans
end

function RadWaveFunc(E,QN::QuantumNumber,A,B,C,rmesh)
    h=rmesh[2]-rmesh[1]

    u=zeros(Float64,Nmesh)
    v=zeros(Float64,Nmesh)
    u[1],v[1],u[Nmesh],v[Nmesh]=BoundCond(QN,E,A[1],B[1],C[1],rmesh)

    for i in 2:Nmatch
        u[i],v[i]=DiffEqBD(u[i-1],v[i-1],E,A[i],B[i],C[i],h)
    end
    u[1:Nmatch]/=u[Nmatch]

    for i in Nmesh-1:-1:Nmatch
        u[i],v[i]=DiffEqBD(u[i+1],v[i+1],E,A[i],B[i],C[i],-h)
    end
    u[Nmatch:Nmesh]/=u[Nmatch]

    Norm=NormFact(rmesh,u)
    u*=sign(u[2]-u[1])/Norm

    return u
end

function CalcStates(QN::QuantumNumber,h2mB,dh2mB,VB,WB,rmesh)
    States=SingleParticleState[]

    A,B,C=CalcABC(QN,h2mB,dh2mB,VB,WB,rmesh)
    Erange=-100.0:0.0 #In BCS approx., E>0 state also needs to be calculated.
    args=[QN,A,B,C,rmesh]
    
    for i in 1:(length(Erange)-1)
        Eans=MyLib.MyBisect(Erange[i],Erange[i+1],WronskyEuler,args,rtol=1e-6) #WronskyEuler(E,QN::QuantumNumber,A,B,C,rmesh)
        if isnan(Eans)==true
            continue
        else
            ψ=RadWaveFunc(Eans,args...)
            push!(States,SingleParticleState(QN,Eans,ψ[:]))
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

function InitialCondition(AN::AtomNum)
    rmesh=getrmesh()
    h2m,dh2m,V,W=InitPot(AN,rmesh)

    InitState=CalcAllStates(h2m,dh2m,V,W,rmesh)
    #Initocc=Calc_occ(AN,InitState)

    return InitState
end

##########################################################
# defene Density, Potential
function Calc_ρ(occ::Vector{Float64},States::Vector{SingleParticleState},rmesh)
    ρ=zeros(Float64,Nmesh)
    for i=eachindex(occ)
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
    dρ=Calc_dρ(ρ,rmesh)
    ddρ=MyLib.diff2nd(h,ρ)
    @. Lapρ[:]+=2*dρ[:]/rmesh[:] + ddρ[:]
    return Lapρ
end

function Calc_τ(occ,States::Vector{SingleParticleState},rmesh)
    τ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    for i=eachindex(occ)

        dRr=MyLib.diff1st(h,(@. States[i].ψ[:]/rmesh[:]))

        j=States[i].QN.j
        l=States[i].QN.l
        @. τ[:]+=occ[i]*(2*j+1)/(4*π)*dRr[:]^2

        if l>0
            @. τ[:]+=occ[i]*(2*j+1)/(4*π)*l*(l+1)*(States[i].ψ[:]/rmesh[:])^2/rmesh[:]^2
        end
    end
    return τ
end

function Calc_J(occ,States::Vector{SingleParticleState},rmesh)
    J=zeros(Float64,Nmesh)
    for i=eachindex(occ)
        j=States[i].QN.j
        l=States[i].QN.l
        if l>0
            @. J[:]+=occ[i]*(2*j+1)/(4*π)*(j*(j+1)-l*(l+1)-0.75)*(States[i].ψ[:]/rmesh[:])^2/rmesh[:]
        end
    end
    return J
end

function Calc_divJ(J::Vector{Float64},rmesh)
    divJ=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    dJ=MyLib.diff1st(h,J)
    @. divJ+=dJ[:]+2*J[:]/rmesh[:]
    return divJ
end

##################################################3
# define potential

function Calc_h2mN(b,aN,aΛ,ρN::Vector{Float64},ρq::Vector{Float64},ρΛ::Vector{Float64})
    QN=QuantumNumber(0,0,b)
    m=getmass(QN)
    return @. ħc^2/(2*m)+aN[5]*ρN[:]+aN[6]*ρq[:]+aΛ[2]*ρΛ[:]
end

function Calc_h2mΛ(aΛ,ρN::Vector{Float64})
    return @. ħc^2/(2*mΛMeV)+aΛ[2]*ρN
end

function Calc_h2mΛ(aΛ,ρN::Vector{Float64})
    return @. ħc^2/(2*mΛMeV)+aΛ[2]*ρN
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
    Vcoul*=e2MeVfm/2 #Chabanat
    #Vcoul*=e2MeVfm #Reainhard

    return Vcoul
    
end

function Calc_Wq(aN,W0,dρN::Vector{Float64},dρq::Vector{Float64},JN::Vector{Float64},Jq::Vector{Float64})
    return @. 0.5*W0*(dρN[:]+dρq[:])+2*aN[9]*JN[:]+2*aN[10]*Jq[:]
end


##########################################################
# HF iteration

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

function CheckConvergence(Oldocc,OldStates,Newocc,NewStates,rmesh;rtol=1e-5)
    check=true
    diffρ=zeros(Float64,3)
    for b in 1:3
        ρold=Calc_ρ(Oldocc[b],OldStates[b],rmesh)
        ρnew=Calc_ρ(Newocc[b],NewStates[b],rmesh)
        #diffρ[b]=MyLib.IntTrap(rmesh,@. (4*π*rmesh[:]^2*(ρold[:]-ρnew[:])))
        diffρ[b]=MyLib.IntTrap(rmesh,(@. 4*π*rmesh[:]^2*(ρold[:]-ρnew[:])^2))
        if abs(diffρ[b])>rtol
            check=false
        end
    end

    println("diffρ=($(diffρ[1]),$(diffρ[2]),$(diffρ[3]))")

    return check

end

function HF_iter(AN::AtomNum;MaxIter=10,NParamType="SLy4",ΛParamType="HPΛ1")
    OldStates=InitialCondition(AN)
    Oldocc=Calc_occ(AN,OldStates)
    rmesh=getrmesh()
    Oldρ3,Olddρ3,OldLapρ3,Oldτ3,OldJ3,OlddivJ3=Calc_Density(Oldocc,OldStates)

    aN=NuclParameters.getaN(NParamType)
    aΛ=LambdaParameters.getaΛ(ΛParamType)
    pN=NuclParameters.getParams(NParamType)
    pΛ=LambdaParameters.getParams(ΛParamType)

    #for debug
    ρptest=zeros(Float64,(MaxIter,Nmesh))

    plot()
    #calc several params
    for i in 1:MaxIter
        #for debug
        ρptest[i,:]=Calc_ρ(Oldocc[1],OldStates[1],rmesh)

        h2m,dh2m,V,W=Calc_Coef(Oldρ3,Oldτ3,OldJ3,aN,aΛ,pN,pΛ,AN.Z)

        NewStates=CalcAllStates(h2m,dh2m,V,W,rmesh)
        Newocc=Calc_occ(AN,NewStates)
        Newρ3,Newdρ3,NewLapρ3,Newτ3,NewJ3,NewdivJ3=Calc_Density(Newocc,NewStates)
        
        if CheckConvergence(Oldocc,OldStates,Newocc,NewStates,rmesh)==true
            return Newocc,NewStates
            break
        end
        println(i)
        OldStates=NewStates
        Oldocc=Newocc
        α=0.5
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


