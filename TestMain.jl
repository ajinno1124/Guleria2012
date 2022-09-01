include("Main.jl")
using Plots

AN=AtomNum(82,126,0) #lead
#AN=AtomNum(8,8,0) #oxigen
#AN=AtomNum(2,2,0) #alpha

function TestInitPot()
    rmesh=getrmesh()
    h2m,dh2m,V,W=InitPot(AN,rmesh)
    plot(xlabel="r",ylabel="V",title="Woods Saxon Z=$(AN.Z), N=$(AN.N)")
    #plot!(rmesh,V[1,:],label="proton")
    #plot!(rmesh,V[2,:],label="neutron")
    #plot!(rmesh,V[3,:],label="Λ")
    plot!(rmesh,(@. W[1,:]/rmesh[:]),label="W/r")
end

function CheckABC()
    rmesh=getrmesh()
    h2m,dh2m,V,W=InitPot(AN,rmesh)

    h2mB=h2m[1,:]
    dh2mB=dh2m[1,:]
    VB=V[1,:]
    WB=W[1,:]

    QN=QuantumNumber(0.5,0,1) #(j,l,B)
    j=QN.j
    l=QN.l

    A=zeros(Float64,Nmesh)
    B=zeros(Float64,Nmesh)
    C=zeros(Float64,Nmesh)
    @. A[:] += -h2mB[:]
    @. B[:] += -dh2mB[:]
    @. C[:] += h2mB[:]*l*(l+1)/(rmesh[:]^2) + VB[:] + dh2mB[:]/rmesh[:] + WB[:]/rmesh[:]*(j*(j+1)-l*(l+1)-0.75)
    
    plot(xlabel="r",ylabel="A,B,C",ylim=(minimum(C)-5,1),title="(j,l,B)=($j,$l,$(QN.B)), Z=$(AN.Z), N=$(AN.N)")
    plot!(rmesh,A,label="A")
    plot!(rmesh,B,label="B")
    plot!(rmesh,C,label="C")
end

function Testgetrmesh(rc)
    h=rc/(Nmesh-0.5)
    rmesh=range(0.5*h,(Nmesh-0.5)*h,length=Nmesh)
    return rmesh
end

function TestInitialCondition()
    
    InitState=InitialCondition(AN)
    Initocc=Calc_occ(AN,InitState)
    for i=eachindex(InitState[1])
        occ=Initocc[1][i]
        l=InitState[1][i].QN.l
        j=InitState[1][i].QN.j
        E=InitState[1][i].E
        println("(occ,l,j,E)=($occ,$l,$j,$E)")
    end
    
    #proton波動関数をプロット
    plot(xlabel="r", ylabel="R", title="Initial Proton Wave Function")
    rmesh=getrmesh()
    for i=eachindex(InitState[1])

            l=InitState[1][i].QN.l
            j=InitState[1][i].QN.j
        if l==0
            plot!(rmesh,InitState[1][i].ψ,label="(l,j)=($l,$j)")
        end
    end
    #println(InitState[1][1].ψ)
    plot!()
end

function TestDensity()
    @time InitState=InitialCondition(AN)
    Initocc=Calc_occ(AN,InitState)
    rmesh=getrmesh()

    @time ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Initocc,InitState)

    Checkρ=MyLib.IntTrap(rmesh,@. 4*π*rmesh[:]^2*ρ3[1,:])
    println("check Z=$(AN.Z): Integrate ρ=$(Checkρ)")
    plot(xlabel="r",ylabel="Density")
    plot!(rmesh,ρ3[1,:],label="ρ")
    plot!(rmesh,dρ3[1,:],label="dρ")
    plot!(rmesh,τ3[1,:],label="τ")
    plot!(rmesh,J3[1,:],label="J")
    plot!(rmesh,Lapρ3[1,:],label="Lapρ")
    plot!(rmesh,divJ3[1,:],label="divJ")
    
end

function CheckCalc_Coef(;NParamType="SLy4",ΛParamType="HPΛ1")
    InitStates=InitialCondition(AN)
    Initocc=Calc_occ(AN,InitStates)
    rmesh=getrmesh()

    aN=NuclParameters.getaN(NParamType)
    aΛ=LambdaParameters.getaΛ(ΛParamType)
    pN=NuclParameters.getParams(NParamType)
    pΛ=LambdaParameters.getParams(ΛParamType)

    h2m,dh2m,V,W=Calc_Coef(Initocc,InitStates,aN,aΛ,pN,pΛ,AN.Z)

    ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Initocc,InitStates)
    ρN=ρ3[1,:]+ρ3[2,:]
    dρN=dρ3[1,:]+dρ3[2,:]
    LapρN=Lapρ3[1,:]+Lapρ3[2,:]
    τN=τ3[1,:]+τ3[2,:]
    JN=J3[1,:]+J3[2,:]
    divJN=divJ3[1,:]+divJ3[2,:]
    h=rmesh[2]-rmesh[1]
    VΛΛ=Calc_VΛΛ(aΛ, pΛ.γ, ρN,LapρN,τN)
    VΛN=Calc_VΛN(aΛ, pΛ.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:])
    VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1,:], τN, τ3[1,:],LapρN,Lapρ3[1,:],divJN,divJ3[1,:])
    VNn=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[2,:], τN, τ3[2,:],LapρN,Lapρ3[2,:],divJN,divJ3[2,:])
    Vcoul=Calc_Vcoul(ρ3[1,:],rmesh,AN.Z)


    for b=1
        if b==1
            h2mconst=ħc^2/(2*mpMeV)
        elseif b==2
            h2mconst=ħc^2/(2*mnMeV)
        elseif b==3
            h2mconst=ħc^2/(2*mΛMeV)
        end
        plot(xlabel="r [fm]", ylabel="[MeV]",title="b=$b $(NParamType) $(ΛParamType) Z=$(AN.Z) N=$(AN.N) Λ=$(AN.Λ)")
        #plot!(rmesh,h2mconst./h2m[b,:],label="h2m")
        #plot!(rmesh,dh2m[b,:],label="dh2m")
        #plot!(rmesh,V[b,:],label="V")
        #plot!(rmesh,VNp,label="VNp")
        #plot!(rmesh,Vcoul,label="Vcoul")
        plot!(rmesh,(@. W[b,:]/rmesh[:]),label="W/r")
    end
    plot!()
end

function CheckABC2(;NParamType="SLy4",ΛParamType="HPΛ1")
    InitStates=InitialCondition(AN)
    Initocc=Calc_occ(AN,InitStates)
    rmesh=getrmesh()

    aN=NuclParameters.getaN(NParamType)
    aΛ=LambdaParameters.getaΛ(ΛParamType)
    pN=NuclParameters.getParams(NParamType)
    pΛ=LambdaParameters.getParams(ΛParamType)

    h2m,dh2m,V,W=Calc_Coef(Initocc,InitStates,aN,aΛ,pN,pΛ,AN.Z)

    h2mB=h2m[1,:]
    dh2mB=dh2m[1,:]
    VB=V[1,:]
    WB=W[1,:]

    QN=QuantumNumber(0.5,0,1) #(j,l,B)
    j=QN.j
    l=QN.l

    A=zeros(Float64,Nmesh)
    B=zeros(Float64,Nmesh)
    C=zeros(Float64,Nmesh)
    @. A[:] += -h2mB[:]
    @. B[:] += -dh2mB[:]
    @. C[:] += h2mB[:]*l*(l+1)/(rmesh[:]^2) + VB[:] + dh2mB[:]/rmesh[:] + WB[:]/rmesh[:]*(j*(j+1)-l*(l+1)-0.75)
    
    plot(xlabel="r",ylabel="A,B,C",title="(j,l,B)=($j,$l,$(QN.B)), Z=$(AN.Z), N=$(AN.N)")
    plot!(rmesh,A,label="A")
    plot!(rmesh,B,label="B")
    plot!(rmesh,C,label="C")
end

function TestHFiter(;NParamType="SLy4",ΛParamType="HPΛ1")
    Ansocc,AnsStates=HF_iter(AN,NParamType=NParamType,ΛParamType=ΛParamType)
    for i=eachindex(AnsStates[1])
        occ=Ansocc[1][i]
        l=AnsStates[1][i].QN.l
        j=AnsStates[1][i].QN.j
        E=AnsStates[1][i].E
        println("(occ,l,j,E)=($occ,$l,$j,$E)")
    end

    rmesh=getrmesh()
    ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Ansocc,AnsStates)
    ρN=ρ3[1,:]+ρ3[2,:]
    dρN=dρ3[1,:]+dρ3[2,:]
    LapρN=Lapρ3[1,:]+Lapρ3[2,:]
    τN=τ3[1,:]+τ3[2,:]
    JN=J3[1,:]+J3[2,:]
    divJN=divJ3[1,:]+divJ3[2,:]
    h=rmesh[2]-rmesh[1]
    aN=NuclParameters.getaN(NParamType)
    aΛ=LambdaParameters.getaΛ(ΛParamType)
    pN=NuclParameters.getParams(NParamType)
    pΛ=LambdaParameters.getParams(ΛParamType)
    VΛΛ=Calc_VΛΛ(aΛ, pΛ.γ, ρN,LapρN,τN)
    VΛN=Calc_VΛN(aΛ, pΛ.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:])
    VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1,:], τN, τ3[1,:],LapρN,Lapρ3[1,:],divJN,divJ3[1,:])
    VNn=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[2,:], τN, τ3[2,:],LapρN,Lapρ3[2,:],divJN,divJ3[2,:])
    Vcoul=Calc_Vcoul(ρ3[1,:],rmesh,AN.Z)
    
    plot(xlabel="r")
    plot!(rmesh,ρ3[1,:],label="ρ")
    plot!(rmesh,VNp,label="VNp")
    plot!(rmesh,Vcoul,label="Vcoul")
    
end