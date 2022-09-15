include("Main.jl")
using Plots


function TestInitPot()
    AN=AtomNum(82,126,0) #lead

    rmesh=getrmesh()
    h2m,dh2m,V,W=InitPot(AN,rmesh)
    labels=["proton" "neutron" "Λ"]
    colors=[:red :blue :green]
    l=@layout [a b c]
    p1=plot(rmesh,[h2m[1,:],h2m[2,:],h2m[3,:]],xlabel="r",ylabel="h2m")
    p2=plot(rmesh,[V[1,:],V[2,:],V[3,:]],xlabel="r",ylabel="V")
    p3=plot(rmesh,[(@. W[1,:]/rmesh[:]),(@. W[2,:]/rmesh[:]),(@. W[3,:]/rmesh[:])],xlabel="r",ylabel="W/r")

    plot(p1,p2,p3,layout=l,title="Z=$(AN.Z), N=$(AN.N)",label=labels,color=colors)
end

function CheckInitABC()
    rmesh=getrmesh()
    AN=AtomNum(82,126,0) #lead
    h2m,dh2m,V,W=InitPot(AN,rmesh)

    A=zeros(Float64,(3,Nmesh))
    B=zeros(Float64,(3,Nmesh))
    C=zeros(Float64,(3,Nmesh))

    for b in 1:3
        QN=QuantumNumber(1.5,1,b)
        A[b,:],B[b,:],C[b,:]=CalcABC(QN,h2m[b,:],dh2m[b,:],V[b,:],W[b,:],rmesh)
    end
    
    labels=["proton" "neutron" "Λ"]
    colors=[:red :blue :green]
    l=@layout [a b c]
    p1=plot(rmesh,[A[1,:],A[2,:],A[3,:]],xlabel="r",ylabel="A")
    p2=plot(rmesh,[B[1,:],B[2,:],B[3,:]],xlabel="r",ylabel="B")
    p3=plot(rmesh,[C[1,:],C[2,:],C[3,:]],xlabel="r",ylabel="C",ylim=(-60,10))

    plot(p1,p2,p3,layout=l,title="Z=$(AN.Z), N=$(AN.N)",label=labels,color=colors)
end

function TestWronskyEuler()
    rmesh=getrmesh()
    AN=AtomNum(82,126,0) #lead
    h2m,dh2m,V,W=InitPot(AN,rmesh)

    A=zeros(Float64,(3,Nmesh))
    B=zeros(Float64,(3,Nmesh))
    C=zeros(Float64,(3,Nmesh))
    b=1 #proton
    QN=QuantumNumber(1.5,1,b)
    A[b,:],B[b,:],C[b,:]=CalcABC(QN,h2m[b,:],dh2m[b,:],V[b,:],W[b,:],rmesh)

    Erange=-100.0:0.1:-15
    Wronskian=zeros(Float64,length(Erange))
    for i=eachindex(Erange)
        Wronskian[i]=WronskyEuler(Erange[i],QN,A[b,:],B[b,:],C[b,:],rmesh)
    end

    plot(xlabel="E (MeV)", ylabel="Wronskian",title="Z=$(AN.Z), N=$(AN.N), j=$(QN.j), l=$(QN.j), b=$(QN.B)")
    plot!(Erange,Wronskian,label=false)

end

function TestInitialCondition()
    AN=AtomNum(82,126,0)
    InitState=InitialCondition(AN)
    Initocc=Calc_occ(AN,InitState)
    b=2
    for i=eachindex(InitState[b])
        occ=Initocc[b][i]
        l=InitState[b][i].QN.l
        j=InitState[b][i].QN.j
        E=InitState[b][i].E
        println("(occ,l,j,E)=($occ,$l,$j,$E)")
    end
    
    #proton波動関数をプロット
    plot(xlabel="r", ylabel="R", title="Initial Proton Wave Function")
    rmesh=getrmesh()
    for i=eachindex(InitState[b])

            l=InitState[b][i].QN.l
            j=InitState[b][i].QN.j
        #if l==0
            #plot!(rmesh,InitState[b][i].ψ,label="(l,j)=($l,$j)",legend=false)
            #plot!(rmesh,(@. InitState[b][i].ψ[:]/rmesh[:]),label="(l,j)=($l,$j)",legend=false)
        #end
        if l==1
            plot!(rmesh,(@. InitState[b][i].ψ[:]^2/rmesh[:]^3),label="(l,j)=($l,$j)",legend=false)
        end
    end
    #println(InitState[1][1].ψ)
    plot!()
end

function TestDensity()
    AN=AtomNum(8,8,0)
    @time InitState=InitialCondition(AN)
    Initocc=Calc_occ(AN,InitState)
    rmesh=getrmesh()

    @time ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Initocc,InitState)

    b=1

    #println(Lapρ3[b,:])
    #println(τ3[b,:])

    Checkρ=MyLib.IntTrap(rmesh,@. 4*π*rmesh[:]^2*ρ3[1,:])
    println("check Z=$(AN.Z): Integrate ρ=$(Checkρ)")
    plot(xlabel="r",ylabel="Density")
    plot!(rmesh,ρ3[b,:],label="ρ")
    plot!(rmesh,dρ3[b,:],label="dρ")
    plot!(rmesh,τ3[b,:],label="τ")
    plot!(rmesh,J3[b,:],label="J")
    plot!(rmesh,Lapρ3[b,:],label="Lapρ")
    #plot!(rmesh,divJ3[b,:],label="divJ")
    
end

function CheckCalc_Coef(;NParamType="SLy4",ΛParamType="HPΛ1")
    AN=AtomNum(2,2,0)
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

    for b=2
        plot(xlabel="r [fm]", ylabel="[MeV]",title="b=$b $(NParamType) $(ΛParamType) Z=$(AN.Z) N=$(AN.N) Λ=$(AN.Λ)")
        plot!(rmesh,h2m[b,:],label="h2m")
        plot!(rmesh,dh2m[b,:],label="dh2m")
        plot!(rmesh,V[b,:],label="V")
        plot!(rmesh,VNp,label="VNp")
        plot!(rmesh,Vcoul,label="Vcoul")
        plot!(rmesh,(@. W[b,:]/rmesh[:]),label="W/r")
    end
    plot!()
end

function TestHFiter(;NParamType="SLy4",ΛParamType="HPΛ1")
    AN=AtomNum(2,2,0)
    Ansocc,AnsStates=HF_iter(AN,NParamType=NParamType,ΛParamType=ΛParamType,MaxIter=50)
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
    plot!()
    
end