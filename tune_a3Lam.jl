include("Main.jl")
using CSV
using DataFrames
using .Threads
const α=0.1

#want to use the same method used in Main.jl
#can be the cause of bugs
function HF_iter_for_a3(AN::AtomNum,aN,aL,pN,pL;MaxIter=50,α=α) #α must be same as run_threads.jl!
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
	rmesh=getrmesh()
	ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Coreocc,CoreStates)
	ρN=ρ3[1,:]+ρ3[2,:]
    dρN=dρ3[1,:]+dρ3[2,:]
    LapρN=Lapρ3[1,:]+Lapρ3[2,:]
    τN=τ3[1,:]+τ3[2,:]
    JN=J3[1,:]+J3[2,:]
    divJN=divJ3[1,:]+divJ3[2,:]

	#=
	E_Kin=Energy_N_Kin(AN,τ3)
	E_SPS=Energy_N_SPS(Coreocc,CoreStates)
	En_R=Energy_N_R(aN,pN.σ,ρ3,ρN)
	Epair=Energy_Pair()
	Etot2=0.5*(E_Kin+E_SPS)- En_R + Epair
	=#
	En=Energy_N(aN,pN.σ,pN.W0,ρ3,ρN,τ3,τN,Lapρ3,LapρN,J3,JN,divJ3,divJN)
	Ec_dir=Energy_coul_dir(ρ3[1,:],rmesh,AN.Z)
	Ec_exch=Energy_coul_exch(ρ3[1,:])
	Epair=Energy_Pair()
	Ecm_dir=Energy_CM_dir(AN,τ3)
	Ecm_exch=Energy_CM_exch()
	Etot = En + Ec_dir + Ec_exch + Epair - Ecm_dir - Ecm_exch

	return Etot

end

function TotalEnergyHYP(Ansocc,AnsStates,AN::AtomNum,aN,aL,pN,pL)
	Etot=0.0
	if  length(AnsStates[3])>=1
		rmesh=getrmesh()
		ρ3,dρ3,Lapρ3,τ3,J3,divJ3=Calc_Density(Ansocc,AnsStates)
		ρN=ρ3[1,:]+ρ3[2,:]
		ρN=ρ3[1,:]+ρ3[2,:]
		dρN=dρ3[1,:]+dρ3[2,:]
		LapρN=Lapρ3[1,:]+Lapρ3[2,:]
		τN=τ3[1,:]+τ3[2,:]
		JN=J3[1,:]+J3[2,:]
		divJN=divJ3[1,:]+divJ3[2,:]
		h=rmesh[2]-rmesh[1]

		#=
		E_N_Kin=Energy_N_Kin(AN,τ3)
		E_L_Kin=Energy_L_Kin(AN,τ3)
		E_N_SPS=Energy_N_SPS(Ansocc,AnsStates)
		En_R=Energy_N_R(aN,pN.σ,ρ3,ρN)
		El_R=Energy_L_R(aL,pL.γ1,pL.γ2,pL.γ3,pL.γ4,ρ3,ρN)
		Epair=Energy_Pair()
		#Etot=0.5*(E_N_Kin+E_N_SPS)- En_R + AnsStates[3][i].E + Epair
		=#
		#Etot=0.5*(E_N_Kin+E_N_SPS)- En_R + Epair + (0.5*(E_L_Kin+AnsStates[3][1].E)-El_R)

		En=Energy_N(aN,pN.σ,pN.W0,ρ3,ρN,τ3,τN,Lapρ3,LapρN,J3,JN,divJ3,divJN)
		Ec_dir=Energy_coul_dir(ρ3[1,:],rmesh,AN.Z)
		Ec_exch=Energy_coul_exch(ρ3[1,:])
		Epair=Energy_Pair()
		Ecm_dir=Energy_CM_dir(AN,τ3)
		Ecm_exch=Energy_CM_exch()
		El=Energy_L(aL,pL.γ1,pL.γ2,pL.γ3,pL.γ4,ρ3,ρN,τ3,τN,Lapρ3,LapρN)
		Etot = En + Ec_dir + Ec_exch + Epair - Ecm_dir - Ecm_exch + El

	else
		Etot+=NaN
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
	Coreocc,CoreStates=HF_iter(AN_core,NParamType=NParamType,LParamType=-1,MaxIter=20,α=α)
	Etot_core=TotalEnergyOfCore(Coreocc,CoreStates,AN_core,aN,pN)

	#println(Etot_core)

	args=(Etot_core,AN,aN,aL,pN,pL,ExpBE)

	# find the crossing Line of a3Lam.
	# calculate 2 times F(a3Lam[i])
	for i in 1:length(a3Lam)-1
		a3Lam_ans=MyLib.MyBisect(a3Lam[i], a3Lam[i+1], EnergyDiff, args;rtol=1e-3)
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
	Etot_HYP=-BE_13LamC+Etot_core

	return a3Lam_ans,BE_13LamC,Etot_HYP,Etot_core

end

function Outputa3()
	df=DataFrame(CSV.File("Lambda Parameters.csv"))
	#index=vcat([32,33,42,43],51:1562)
	#index=1:3578
	index=47:50

	AN=AtomNum(6,6,1)
	ExpBE=11.88 #MeV, Hashimoto & Tamura, Gogami
	NParamType="SLy4"
	a3Lam_ans=zeros(Float64,length(index))
	BE_13LamC=zeros(Float64,length(index))
	Etot_HYP=zeros(Float64,length(index))
	Etot_core=zeros(Float64,length(index))
	@threads for i=eachindex(index)
		println("index = $(index[i])")
		a3Lam_ans[i],BE_13LamC[i],Etot_HYP[i],Etot_core[i]=tune_a3Lam(AN,ExpBE,NParamType=NParamType,LParamType=index[i])
		println("index = $(index[i]), a3Lam_ans=$(a3Lam_ans[i]), BE_13LamC=$(BE_13LamC[i]) (MeV)")
	end

	io1=open("a3.csv","w")
	write(io1,"index,Parameter Name,a3,BE(13LamC)(MeV),BE-11.88MeV,Etot_HYP(MeV),Etot_core(MeV)\n")
	for i=eachindex(index)
		write(io1,"$(df[index[i],"index"])")
		write(io1,",$(df[index[i],"Parameter Name"])")
		write(io1,",$(a3Lam_ans[i])")
		write(io1,",$(BE_13LamC[i])") #calculate again
		write(io1,",$(BE_13LamC[i]-11.88)")
		write(io1,",$(Etot_HYP[i])")
		write(io1,",$(Etot_core[i])\n")
	end
	close(io1)
end

@time Outputa3()