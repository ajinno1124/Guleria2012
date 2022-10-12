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
	return (@. p[1]+p[2]*k[:]^2)
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

function VLL(ρp,ρn,aL1,aL2,aL4,aL5)
	ρN=ρp+ρn
    τN=getτq.(ρp)+getτq.(ρn)

	return (@. aL1*ρN + aL2*τN + aL4*ρN^(4/3) + aL4*ρN^(5/3))
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
    τN=getτq.(ρp)+getτq.(ρn)

	aL1=C
	aL1/=ρN - aL2*τN - aL4*ρN^(4/3) - aL5*ρN^(5/3)
	return aL1
end

function get_aL1_PNM(C,aL2,aL4,aL5)
	ρ0=ρ0_Kohno
	ρp=0.0
	ρn=ρ0
    τN=getτq.(ρp)+getτq.(ρn)

	aL1=C
	aL1/=ρN - aL2*τN - aL4*ρN^(4/3) - aL5*ρN^(5/3)
	return aL1
end

function getaL2(Paramtype::String)
	C,aL2=[0.0,0.0]
	if ParamType=="MD1"
	elseif ParamType=="MD2_2.5"
		C,aL2=[-28.278161390790633, 5.97294433002585]
	elseif ParamType=="MD3_1.0"
		C,aL2=[-29.95328507758762, 6.79868183475533]
	end

	return C,aL2

end

function Fit_Dense(dataname;isSNM=false,isPNM=false,ParamType)
	C,aL2=getaL2(Paramtype)
	data,header=readdlm(dataname,',',header=true)
    df=DataFrame(data, vec(header))
    xdata=ρ0*df[df.density.<ρ_cutoff,:].density
    ydata=df[df.density.<ρ_cutoff,:].U

    p0=[0.0,0.0,0.0]
    if isSNM==true
		get_aL1_SNM(C,aL2,aL4,aL5)
		Func(ρN,p)=FitVΛΛ_SNM(ρN,,,aL4,aL5)
        fit=curve_fit(FitVΛΛ_SNM, xdata, ydata, p0)
    elseif isPNM==true
		Func(ρN,)
        fit=curve_fit(FitVΛΛ_PNM, xdata, ydata, p0)
    end

    println("\ndegree of freedom = $(dof(fit))") #degrees of freedom
    println("aΛ1,aΛ2,aΛ4 = $(coef(fit))") #best fit parameters
    #println("residual = $(fit.resid)") #residuals = vector of residuals
	#println("jacobian = $(fit.jacobian)") #estimated Jacobian at solution

    plot(xlabel="ρ/ρ0", ylabel="U",title="$(dataname)")
    scatter!(xdata/ρ0, ydata, label="data")
    plot!(df.density, df.U, label="all data")
    ρmesh=0:0.1:5
    if isSNM==true
        plot!(ρmesh,FitVΛΛ_SNM(ρmesh*ρ0,coef(fit)),label="fit")
    elseif isPNM==true
        plot!(ρmesh,FitVΛΛ_PNM(ρmesh*ρ0,coef(fit)),label="fit")
    end
    plot!()
end