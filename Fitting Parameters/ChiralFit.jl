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
    ρ0=0.16

    ρ_cutoff=1.0
    k_cutoff=2.0
end


function getkfq(ρq)
    return (6*π^2*ρq/gq)^(1/3)
end

function getτq(ρq)
    kfq=getkfq(ρq)
    return 3.0/5.0*kfq^2*ρq
end

function VΛΛ(ρp,ρn,aΛ1,aΛ2,aΛ4)
    γ=1/3
    ρN=ρp+ρn
    τN=getτq.(ρp)+getτq.(ρn)

    return @. aΛ1*ρN+aΛ2*τN+aΛ4*ρN^(γ+1)
end

function FitVΛΛ_SNM(ρN,params::Vector{Float64})
    return VΛΛ(ρN/2.0,ρN/2.0,params...)
end

function FitVΛΛ_PNM(ρN,params::Vector{Float64})
    return VΛΛ(0.0,ρN,params...)
end

function FittingMain(dataname;isSNM=false,isPNM=false)
    data,header=readdlm(dataname,',',header=true)
    df=DataFrame(data, vec(header))
    xdata=ρ0*df[df.density.<ρ_cutoff,:].density
    ydata=df[df.density.<ρ_cutoff,:].U

    p0=[0.0,0.0,0.0]
    if isSNM==true
        fit=curve_fit(FitVΛΛ_SNM, xdata, ydata, p0)
    elseif isPNM==true
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

###########################################3
# Calc LParams, from ,a, b, c

# aΛ[3]=0.125*(3u1-u2)=0
function Calc_LParams1(a,b,c)
    γ=1/3
    u0=a/ρ0
    u1=c/ρ0^(5/3)*5.0/3.0*(3*π^2/2)^(-2/3)
    u2=3*u1
    u3=0
    u3p=8.0/3.0*b/ρ0^(4/3)
    y0=0
    y3=0
    return [γ,u0,u1,u2,u3,u3p,y0,y3]
end

function Calc_GKWCoef1()
    a,b,c=[-154.9, 142.4, -21.4]
    println("GKW2: [γ,u0,u1,u2,u3,u3p,y0,y3]")
    println("= $(Calc_LParams1(a,b,c))\n")
    a,b,c=[-80.1, 0.16, 50.4]
    println("GKW3: [γ,u0,u1,u2,u3,u3p,y0,y3]")
    println("= $(Calc_LParams1(a,b,c))\n")
end

include("../SkyrmeParams.jl")
using ..LambdaParameters
function ComparePot()
    labels=["HPL2","HPL3","NL1","OL1","GKW2", "GKW3"]
    ρmesh=0:0.1:5
    plot(xlabel="ρ/ρ0", ylabel="UΛ",legend=:topleft)
    for i=eachindex(labels)
        aΛ=LambdaParameters.getaΛ(labels[i])
        UΛ=FitVΛΛ_SNM(ρmesh*ρ0,[aΛ[1],aΛ[2],aΛ[4]])
        plot!(ρmesh,UΛ,label=labels[i])
    end
    plot!()

end