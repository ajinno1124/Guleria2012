include("Main.jl")
using Plots


function TestInitPot()
    AN=AtomNum(82,126,0) #lead

    rmesh=getrmesh()
    h2m,dh2m,ddh2m,V,W=InitPot(AN,rmesh)
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
    h2m,dh2m,ddh2m,V,W=InitPot(AN,rmesh)

    A=zeros(Float64,(3,Nmesh))
    C=zeros(Float64,(3,Nmesh))

    for b in 1:3
        QN=QuantumNumber(1.5,1,b)
        A[b,:],C[b,:]=CalcABC(QN,h2m[b,:],dh2m[b,:],ddh2m[b,:],V[b,:],W[b,:],rmesh)
    end
    
    labels=["proton" "neutron" "Λ"]
    colors=[:red :blue :green]
    l=@layout [a c]
    p1=plot(rmesh,[A[1,:],A[2,:],A[3,:]],xlabel="r",ylabel="A")
    p2=plot(rmesh,[C[1,:],C[2,:],C[3,:]],xlabel="r",ylabel="C",ylim=(-60,10))

    plot(p1,p2,layout=l,title="Z=$(AN.Z), N=$(AN.N)",label=labels,color=colors)
end

function TestWronskyEuler()
    rmesh=getrmesh()
    AN=AtomNum(82,126,0) #lead
    h2m,dh2m,ddh2m,V,W=InitPot(AN,rmesh)

    A=zeros(Float64,(3,Nmesh))
    C=zeros(Float64,(3,Nmesh))
    b=1 #proton
    QN=QuantumNumber(0.5,0,b)
    A[b,:],C[b,:]=CalcABC(QN,h2m[b,:],dh2m[b,:],ddh2m[b,:],V[b,:],W[b,:],rmesh)

    Erange=-65.0:0.1:-15
    Wronskian=zeros(Float64,length(Erange))
    for i=eachindex(Erange)
        Wronskian[i]=WronskyEuler(Erange[i],QN,A[b,:],C[b,:],rmesh)
    end

    plot(xlabel="E (MeV)", ylabel="Wronskian",title="Z=$(AN.Z), N=$(AN.N), j=$(QN.j), l=$(QN.l), b=$(QN.B)")
    plot!(Erange,Wronskian,label=false)

end

function TestInitialCondition()
    AN=AtomNum(82,126,0)
    InitState=InitialCondition(AN)
    Initocc=Calc_occ(AN,InitState)
    b=1
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
            plot!(rmesh,(@. InitState[b][i].ψ[:]),label="(l,j)=($l,$j)",legend=false)
        #end
        #if l==1
        #    plot!(rmesh,(@. InitState[b][i].ψ[:]^2/rmesh[:]^2),label="(l,j)=($l,$j)",legend=false)
        #end
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
    #println(dρ3[b,:])
    #println(ρ3[b,:])

    Checkρ=MyLib.IntTrap(rmesh,@. 4*π*rmesh[:]^2*ρ3[1,:])
    Checkρ+=4*π*rmesh[1]^3*ρ3[1,1]/2
    println("check Z=$(AN.Z): Integrate ρ=$(Checkρ)")
    plot(xlabel="r",ylabel="Density",xlim=(0,10))
    plot!(rmesh,ρ3[b,:],label="ρ")
    plot!(rmesh,dρ3[b,:],label="dρ")
    plot!(rmesh,τ3[b,:],label="τ")
    plot!(rmesh,J3[b,:],label="J")
    plot!(rmesh,Lapρ3[b,:],label="Lapρ")
    plot!(rmesh,divJ3[b,:],label="divJ")

end

function TestHFiter(;AN=AtomNum(82,125,1),NParamType="SLy4",LParamType="HPL2",α=0.5)
    #AN=AtomNum(82,125,1)
    #AN=AtomNum(82,126,1)
    Ansocc,AnsStates=HF_iter(AN,NParamType=NParamType,LParamType=LParamType,MaxIter=50,α=α)
    b=3
    for i=eachindex(AnsStates[b])
        occ=Ansocc[b][i]
        l=AnsStates[b][i].QN.l
        j=AnsStates[b][i].QN.j
        E=AnsStates[b][i].E
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

	dτ3=zeros(Float64,(3,Nmesh))
    ddρ3=zeros(Float64,(3,Nmesh))
    for b in 1:3
        dτ3[b,:]=Calc_dτ(τ3[b,:],rmesh)
        ddρ3[b,:]=Calc_ddρ(ρ3[b,:],rmesh)
    end
    dτN=dτ3[1,:]+dτ3[2,:]
    ddρN=ddρ3[1,:]+ddρ3[2,:]

    aN=NuclParameters.getaN(NParamType)
    aL=LambdaParameters.getaL(LParamType)
    pN=NuclParameters.getParams(NParamType)
    pL=LambdaParameters.getParams(LParamType)
	# Guleria
	if isGuleria==1
		VΛΛ=Calc_VΛΛ_G(aL,pL.γ1,ρN,ddρN,LapρN,τN,dτ3[3,:])
		VΛp=Calc_VΛN_G(aL,pL.γ1,ρ3[3,:],ddρ3[3,:],τ3[3,:],dτN,Lapρ3[3,:],ρN)
		VΛn=VΛp
	else
		# Rayet
		VΛΛ=Calc_VΛΛ(aL, pL.γ1, pL.γ2, pL.γ3, pL.γ4, ρN,LapρN,τN,ρ3[1,:],ρ3[2,:])
		VΛp=Calc_VΛN(aL, pL.γ1, pL.γ2, pL.γ3, pL.γ4, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:],ρ3[1,:])
		VΛn=Calc_VΛN(aL, pL.γ1, pL.γ2, pL.γ3, pL.γ4, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:],ρ3[2,:])
	end
    #VΛΛ=Calc_VΛΛ(aL, pL.γ, ρN,LapρN,τN,ρ3[1,:],ρ3[2,:])
    #VΛp=Calc_VΛN(aL, pL.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:],ρ3[1,:])
    #VΛn=Calc_VΛN(aL, pL.γ, ρN, ρ3[3,:],Lapρ3[3,:],τ3[3,:],ρ3[2,:])
    VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1,:], τN, τ3[1,:],LapρN,Lapρ3[1,:],divJN,divJ3[1,:])
    VNn=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[2,:], τN, τ3[2,:],LapρN,Lapρ3[2,:],divJN,divJ3[2,:])
    Vcoul=Calc_Vcoul(ρ3[1,:],rmesh,AN.Z)

    plot(xlabel="r")
    #plot!(rmesh,ρ3[1,:],label="proton")
    #plot!(rmesh,ρ3[2,:],label="neutron")
    #plot!(rmesh,ρ3[3,:],label="Λ")
    plot!(rmesh,VNp,label="VNp")
    plot!(rmesh,VΛΛ,label="VΛΛ",xlim=(0,10))
    #plot!(rmesh,VΛN,label="VΛN",xlim=(0,10))
    #plot!(rmesh,Vcoul,label="Vcoul")
    plot!()

end

function Test138La()
    AN=AtomNum(57,81,0)
    #α=[0.01,0.05,0.1,0.2,0.3,0.4,0.5] 0.2~0.5 doesn't converge
    α=0.01:0.01:0.15

    NParamType="SLy4"
    LParamType=-1
    #LParamType=1


    aN=NuclParameters.getaN(NParamType)
    aL=LambdaParameters.getaL(LParamType)
    pN=NuclParameters.getParams(NParamType)
    pL=LambdaParameters.getParams(LParamType)

    plot(xlabel="r (fm)", ylabel="VNp",xlim=(0,10),ylim=(-80,0))
    for i=eachindex(α)
        Ansocc,AnsStates,Check=HF_iter(AN,NParamType=NParamType,LParamType=LParamType,MaxIter=50,α=α[i])
        if Check==false
            println("Using α=$(α[i]), do not converge.")
        else
            rmesh=getrmesh()
            ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Ansocc,AnsStates)
            ρN=ρ3[1,:]+ρ3[2,:]
            dρN=dρ3[1,:]+dρ3[2,:]
            LapρN=Lapρ3[1,:]+Lapρ3[2,:]
            τN=τ3[1,:]+τ3[2,:]
            JN=J3[1,:]+J3[2,:]
            divJN=divJ3[1,:]+divJ3[2,:]
            h=rmesh[2]-rmesh[1]

            VNp=Calc_VNq(aN, pN.σ, pN.W0, ρN, ρ3[1,:], τN, τ3[1,:],LapρN,Lapρ3[1,:],divJN,divJ3[1,:])

            plot!(rmesh,VNp,label="α=$(α[i])")
        end
    end

    plot!()

end