#fit the optcal potential
# U_opt = aL1*ρN + aL2*(k^2*ρN+τN) + aL4*ρN^(4/3) + aL5*ρN^(5/3)


using Parameters
using LsqFit
using Plots
using DataFrames
using DelimitedFiles

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

    ρ_cutoff=3.0
	k_cutoff=2.5
end

# fit Mom dependence by U = const + a2*k^2
function Func_fitmom(k,p::Vector{Float64})
	ρ0=ρ0_Kohno
	return (@. p[1]+p[2]*k[:]^2*ρ0)
end

function Fit_Mom(dataname)
	println("k_cutoff=$(k_cutoff)")

    data,header=readdlm(dataname,',',header=true)
    df=DataFrame(data, vec(header))
    xdata=df[df.k.<k_cutoff,:].k
    ydata=df[df.k.<k_cutoff,:].Um

    p0=[0.0,0.0]
	fit=curve_fit(Func_fitmom, xdata, ydata, p0)

    println("\ndegree of freedom = $(dof(fit))") #degrees of freedom
    println("const, aΛ2 = $(coef(fit))") #best fit parameters
    #println("residual = $(fit.resid)") #residuals = vector of residuals
	#println("jacobian = $(fit.jacobian)") #estimated Jacobian at solution

    plot(xlabel="k", ylabel="Um",title="$(dataname)")
    scatter!(df.k, df.Um, label="data")
    kmesh=0:0.1:3
	plot!(kmesh,Func_fitmom(kmesh,coef(fit)),label="fit")
    plot!()

end

function getkfq(ρq)
    return (6*π^2*ρq/gq)^(1/3)
end

function getτq(ρq)
    kfq=getkfq(ρq)
    return 3.0/5.0*kfq^2*ρq
end

function VΛΛ(ρp,ρn,aL1,aL2,aL4,aL5)
	ρN=zeros(Float64,length(ρn))
	@. ρN=ρp+ρn
    τN=getτq.(ρp)+getτq.(ρn)

	return (@. aL1*ρN + aL2*τN + aL4*ρN^(4/3) + aL5*ρN^(5/3))
end

function FitVΛΛ_SNM(ρN,aL1,aL2,aL4,aL5)
    return VΛΛ(ρN/2.0,ρN/2.0,aL1,aL2,aL4,aL5)
end

function FitVΛΛ_PNM(ρN,aL1,aL2,aL4,aL5)
    return VΛΛ(0.0,ρN,aL1,aL2,aL4,aL5)
end

function get_aL1_SNM(C,aL2,aL4,aL5)
	ρ0=ρ0_Kohno
	ρp=ρ0/2.0
	ρn=ρp
	ρN=ρp+ρn
    τN=getτq.(ρp)+getτq.(ρn)

	aL1=C
	aL1/=ρN - aL2*τN - aL4*ρN^(4/3) - aL5*ρN^(5/3)
	return aL1
end

function get_aL1_PNM(C,aL2,aL4,aL5)
	ρ0=ρ0_Kohno
	ρp=0.0
	ρn=ρ0
	ρN=ρp+ρn
    τN=getτq(ρp)+getτq(ρn)

	aL1=C
	aL1/=ρN - aL2*τN - aL4*ρN^(4/3) - aL5*ρN^(5/3)
	return aL1
end

function getaL2(ParamType::String)
	C,aL2=[0.0,0.0]
	if ParamType=="MD1_2.5"
		C,aL2=[-29.098371987745747, 32.293882105765896]
	elseif ParamType=="MD2_2.5"
		C,aL2=[-28.278161390800726, 35.9815923495766]
	elseif ParamType=="MD3_1.0"
		C,aL2=[-29.95328507758867, 40.95591466722299]
	end

	return C,aL2

end

function Fit_Dense(dataname;isSNM=false,isPNM=false,ParamType::String)
	ρ0=ρ0_GKW
	C,aL2=getaL2(ParamType)
	data,header=readdlm(dataname,',',header=true)
    df=DataFrame(data, vec(header))
    xdata=ρ0*df[df.density.<ρ_cutoff,:].density
    ydata=df[df.density.<ρ_cutoff,:].U

    p0=[0.0,0.0]
    if isSNM==true
		Func(ρN::Vector{Float64},p)=FitVΛΛ_SNM(ρN,get_aL1_SNM(C,aL2,p[1],p[2]),aL2,p[1],p[2])
        fit=curve_fit(Func, xdata, ydata, p0)
    elseif isPNM==true
		Func(ρN,p)=FitVΛΛ_PNM(ρN,get_aL1_PNM(C,aL2,p[1],p[2]),aL2,p[1],p[2])
        fit=curve_fit(Func, xdata, ydata, p0)
    end

    println("\ndegree of freedom = $(dof(fit))") #degrees of freedom
    println("aΛ4, aΛ5 = $(coef(fit))") #best fit parameters
    #println("residual = $(fit.resid)") #residuals = vector of residuals
	#println("jacobian = $(fit.jacobian)") #estimated Jacobian at solution

    plot(xlabel="ρ/ρ0", ylabel="U",title="$(dataname)")
    scatter!(xdata/ρ0, ydata, label="used data")
    plot!(df.density, df.U, label="all data")
    umesh=0:0.1:5
    if isSNM==true
		aL1=get_aL1_SNM(C,aL2,coef(fit)[1],coef(fit)[2])
		println("aΛ1, aΛ2, aΛ4, aΛ5 = [$(aL1), $(aL2), $(coef(fit)[1]), $(coef(fit)[2])]")
        plot!(umesh,FitVΛΛ_SNM(umesh*ρ0,aL1,aL2,coef(fit)[1],coef(fit)[2]),label="fit")
    elseif isPNM==true
        aL1=get_aL1_PNM(C,aL2,coef(fit)[1],coef(fit)[2])
		println("aΛ1, aΛ2, aΛ4, aΛ5 = [$(aL1), $(aL2), $(coef(fit)[1]), $(coef(fit)[2])]")
        plot!(umesh,FitVΛΛ_PNM(umesh*̢ρ0,aL1,aL2,coef(fit)[1],coef(fit)[2]),label="fit")
    end
    plot!()
end