include("Main.jl")
include("SkyrmeParams.jl")
@consts begin
    gq=2 #spin degeneracy
    # Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
    mpMeV=938.27208816
    mnMeV=939.5654205
    mAveMeV=(mpMeV+mnMeV)/2
    mΛMeV=1115.683 #\pm 0.006 MeV
    ħc=197.3269804
    ρ0_Kohno=0.166
	ρ0_GKW=0.16
end

function getkfq(ρq)
    return (6*π^2*ρq/gq)^(1/3)
end

function getτq(ρq)
    kfq=getkfq(ρq)
    return 3.0/5.0*kfq^2*ρq
end

function Uopt(ρp,ρn,kL,aL,γ1,γ2)
	ρN=ρp+ρn
	τN=getτq(ρp)+getτq(ρn)
	return aL[1]*ρN+aL[2]*(kL^2*ρN+τN)+aL[4]*ρN^(1+γ1)+aL[5]*ρN^(1+γ2)+aL[6]*(ρN^2 + 2*ρp*ρn)
end

function OutputPotentialData(LParamType::String)
	ρ0=0.16
	pL=LambdaParameters.getParams(LParamType)
	aL=LambdaParameters.getaL(LParamType)

	umesh=0:0.01:5
	#U_dense=zeros(Float64,length(umesh))
	#for i=eachindex(umesh)
	#	U_dense[i]=Uopt(umesh[i]*ρ0/2.0,umesh[i]*ρ0/2.0,0.0,aL,pL.γ1,pL.γ2)
	#end
	Uopt_dense(ρN)=Uopt(ρN/2.0,ρN/2.0,0.0,aL,pL.γ1,pL.γ2)
	U_dense=Uopt_dense.(umesh*ρ0)

	kmesh=0:0.01:3
	Uopt_mom(kL)=Uopt(ρ0/2.0,ρ0/2.0,kL,aL,pL.γ1,pL.γ2)
	U_mom=Uopt_mom.(kmesh)

	rm("data/Potential_$(LParamType)",force=true,recursive=true)
    mkpath("data/Potential_$(LParamType)")
    cd("data/Potential_$(LParamType)")

	io1=open("DensityDep_$(LParamType).csv","w")
	write(io1,"density,U\n")
	for i=eachindex(umesh)
		write(io1,"$(umesh[i]),$(U_dense[i])\n")
	end
	close(io1)

	io2=open("MomentumDep_$(LParamType).csv","w")
	write(io2,"k,Um\n")
	for i=eachindex(kmesh)
		write(io2,"$(kmesh[i]),$(U_mom[i])\n")
	end
	close(io2)


	cd("../..")
end

OutputPotentialData("HPL2")
OutputPotentialData("LY1")
#OutputPotentialData("GKW2")
#OutputPotentialData("GKW3")
#OutputPotentialData("GKW2+MD1")
#OutputPotentialData("GKW3+MD2")
#OutputPotentialData("GKW3+MD3")
OutputPotentialData("GKW2_1.5")
OutputPotentialData("GKW3_1.5")
OutputPotentialData("GKW2_1.5+Kohno2")
OutputPotentialData("GKW3_1.5+Kohno3")