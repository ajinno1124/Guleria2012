include("Main.jl")
using CSV
using DataFrames
using .Threads

#want to use the same method used in Main.jl
#can be the cause of bugs
function HF_iter_for_a3(AN::AtomNum,aN,aL,pN,pL;MaxIter=15,α=0.5)
    OldStates=InitialCondition(AN)
    Oldocc=Calc_occ(AN,OldStates)
    rmesh=getrmesh()
    Oldρ3,Olddρ3,OldLapρ3,Oldτ3,OldJ3,OlddivJ3=Calc_Density(Oldocc,OldStates)

    #aN=NuclParameters.getaN(NParamType)
    #aL=LambdaParameters.getaL(LParamType)
    #pN=NuclParameters.getParams(NParamType)
    #pL=LambdaParameters.getParams(LParamType)

    #for debug
    #ρptest=zeros(Float64,(MaxIter,Nmesh))

    #plot()
    #calc several params
    for i in 1:MaxIter
        #for debug
        #ρptest[i,:]=Calc_ρ(Oldocc[1],OldStates[1],rmesh)

        h2m,dh2m,ddh2m,V,W=Calc_Coef(Oldρ3,Oldτ3,OldJ3,aN,aL,pN,pL,AN)

        NewStates=CalcAllStates(h2m,dh2m,ddh2m,V,W,rmesh)
        Newocc=Calc_occ(AN,NewStates)
        Newρ3,Newdρ3,NewLapρ3,Newτ3,NewJ3,NewdivJ3=Calc_Density(Newocc,NewStates)

        if CheckConvergence(Oldocc,OldStates,Newocc,NewStates,rmesh)==true
            return Newocc,NewStates
            break
        end
        println(i)
        OldStates=NewStates
        Oldocc=Newocc
        Oldρ3=Oldρ3*(1-α)+Newρ3*α
        Oldτ3=Oldτ3*(1-α)+Newτ3*α
        OldJ3=OldJ3*(1-α)+NewJ3*α
    end

    #=
    plot(xlabel="r",ylabel="ρn")
    for i in 1:MaxIter
        plot!(rmesh,ρptest[i,:],label="$i")
    end
    plot!()
    =#

end

function TotalEnergyOfCore(Coreocc,CoreStates,AN::AtomNum,aN,pN)
	ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Coreocc,CoreStates)
	ρN=ρ3[1,:]+ρ3[2,:]

	E_Kin=Energy_N_Kin(AN,τ3)
	E_SPS=Energy_N_SPS(Coreocc,CoreStates)
	En_R=Energy_N_R(aN,pN.σ,ρ3,ρN)
	Epair=Energy_Pair()
	Etot2=0.5*(E_Kin+E_SPS)- En_R + Epair

	return Etot2

end

function TotalEnergyHYP(Ansocc,AnsStates,AN::AtomNum,aN,aL,pN,pL)
    ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Ansocc,AnsStates)
    ρN=ρ3[1,:]+ρ3[2,:]

	E_N_Kin=Energy_N_Kin(AN,τ3)
	E_L_Kin=Energy_L_Kin(AN,τ3)
	E_N_SPS=Energy_N_SPS(Ansocc,AnsStates)
	En_R=Energy_N_R(aN,pN.σ,ρ3,ρN)
	El_R=Energy_L_R(aL,pL.γ1,pL.γ2,pL.γ3,pL.γ4,ρ3,ρN)
	Epair=Energy_Pair()
	#Etot=0.5*(E_N_Kin+E_N_SPS)- En_R + AnsStates[3][i].E + Epair
	if  length(AnsStates[3])>=1
		Etot=0.5*(E_N_Kin+E_N_SPS)- En_R + Epair + (0.5*(E_L_Kin+AnsStates[3][1].E)-El_R)
	else
		Etot=NaN
	end

	return Etot
end

# Calculate Binding Energy for given AN, NParamType, and LParamType
# BENaN=Total Energy of Core Nucleus
function LamBindingEnergy(a3Lam,Etot_core,AN::AtomNum,aN,aL,pN,pL)
	aL_changed=aL
	pL_changed=pL
	aL_changed[3]=a3Lam
	pL_changed.u1=aL[2]+2*aL[3]
	pL_changed.u2=3*aL[2]-2*aL[3]

	Ansocc,AnsStates=HF_iter_for_a3(AN,aN,aL_changed,pN,pL_changed)

	EtotHYP=TotalEnergyHYP(Ansocc,AnsStates,AN,aN,aL_changed,pN,pL_changed)

	return Etot_core-EtotHYP
end

function EnergyDiff(a3Lam,Etot_core,AN::AtomNum,aN,aL,pN,pL,GivenBE)
	BELam=LamBindingEnergy(a3Lam,Etot_core,AN,aN,aL,pN,pL)

	#println("Binding Energy of Lambda for a3Lam=$(a3Lam) is $(BELam) MeV")

	return BELam-GivenBE
end

function tune_a3Lam(AN::AtomNum,ExpBE;NParamType,LParamType)
	#a3Lam=-10:10:200 #possible value of a3Lam
	a3Lam=-10:50:190
	a3Lam_ans=NaN

	aN=NuclParameters.getaN(NParamType)
    aL=LambdaParameters.getaL(LParamType)
    pN=NuclParameters.getParams(NParamType)
    pL=LambdaParameters.getParams(LParamType)

	AN_core=AtomNum(AN.Z,AN.N,0)
	Coreocc,CoreStates=HF_iter(AN_core,NParamType=NParamType,LParamType=-1,MaxIter=20)
	Etot_core=TotalEnergyOfCore(Coreocc,CoreStates,AN_core,aN,pN)

	#println(Etot_core)

	args=(Etot_core,AN,aN,aL,pN,pL,ExpBE)

	# find the crossing Line of a3Lam.
	# calculate 2 times F(a3Lam[i])
	for i in 1:length(a3Lam)-1
		a3Lam_ans=MyLib.MyBisect(a3Lam[i], a3Lam[i+1], EnergyDiff, args;rtol=1e-5)
		if isnan(a3Lam_ans)==true
            continue
        else
            break
        end
	end

	if isnan(a3Lam_ans)
		println("No solution for a3Lam is found.")
	end

	#BE_13LamC=EnergyDiff(a3Lam_ans,Etot_core,AN,aN,aL,pN,pL,GivenBE)
	BE_13LamC=LamBindingEnergy(a3Lam_ans,Etot_core,AN,aN,aL,pN,pL)

	return a3Lam_ans,BE_13LamC

end

#AN=AtomNum(6,6,1)
#@time tune_a3Lam(AN,11.88,NParamType="SLy4",LParamType=33)

function Outputa3()
	df=DataFrame(CSV.File("Lambda Parameters.csv"))
	#index=vcat([32,33,42,43],51:1562)
	#index=vcat([32,33,42,43],51:300)
	#index=240:248
	#index=301:1200
	#index=1201:1562
	index=1563:3578
	#index=1

	AN=AtomNum(6,6,1)
	ExpBE=11.88 #MeV, Hashimoto & Tamura, Gogami
	NParamType="SLy4"
	a3Lam_ans=zeros(Float64,length(index))
	BE_13LamC=zeros(Float64,length(index))
	@threads for i=eachindex(index)
		println("index = $(index[i])")
		a3Lam_ans[i],BE_13LamC[i]=tune_a3Lam(AN,ExpBE,NParamType=NParamType,LParamType=index[i])
		println("index = $(index[i]), a3Lam_ans=$(a3Lam_ans[i]), BE_13LamC=$(BE_13LamC[i])")
	end

	io1=open("a3.csv","w")
	write(io1,"index,Parameter Name,a3,BE(13LamC)(MeV),BE-11.88MeV\n")
	for i=eachindex(index)
		write(io1,"$(df[index[i],"index"])")
		write(io1,",$(df[index[i],"Parameter Name"])")
		write(io1,",$(a3Lam_ans[i])")
		write(io1,",$(BE_13LamC[i])") #calculate again
		write(io1,",$(BE_13LamC[i]-11.88)\n")
	end
	close(io1)
end

@time Outputa3()