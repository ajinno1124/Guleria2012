include("SkyrmeParams.jl")
include("MyLib.jl")
using Parameters
using CSV
using DataFrames
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

    NumberOfParams=50
end

function getkfq(ρq)
    return (6*π^2*ρq/gq)^(1/3)
end

function getτq(ρq)
    kfq=getkfq(ρq)
    return 3.0/5.0*kfq^2*ρq
end

function Uopt(ρp,ρn,kL,aL,γ1,γ2,γ3,γ4)
	ρN=ρp+ρn
	τN=getτq(ρp)+getτq(ρn)
	return aL[1]*ρN+aL[2]*(kL^2*ρN+τN)+aL[4]*ρN^(1+γ1)+aL[5]*ρN^(1+γ2)+aL[6]*ρN^(1+γ3)+aL[7]*ρN^(1+γ4)+aL[8]*(ρN^2+2*ρp*ρn)
end

function Calc_J(aL,γ1,γ2,γ3,γ4)
    ρ0=ρ0_GKW
    τp=getτq(ρ0*0.5)
    τN=2*τp
    J=aL[1]*ρ0+aL[2]*τN
    J+=aL[4]*ρ0^(γ1+1)+aL[5]*ρ0^(γ2+1)+aL[6]*ρ0^(γ3+1)+aL[7]*ρ0^(γ4+1)
    J+=aL[8]*3/2*ρ0^2
    return J
end

function Calc_L(aL,γ1,γ2,γ3,γ4)
    ρ0=ρ0_GKW
    τp=getτq(ρ0*0.5)
    τN=2*τp
    L=aL[1]*3*ρ0+aL[2]*5*τN
    L+=3*( aL[4]*(γ1+1)*ρ0^(γ1+1)+aL[5]*(γ2+1)*ρ0^(γ2+1)+aL[6]*(γ3+1)*ρ0^(γ3+1)+aL[7]*(γ4+1)*ρ0^(γ4+1) )
    L+=aL[8]*9*ρ0^2
    return L
end

function Calc_K(aL,γ1,γ2,γ3,γ4)
    ρ0=ρ0_GKW
    τp=getτq(ρ0*0.5)
    τN=2*τp
    K=aL[2]*10*τN
    K+=9*( aL[4]*(γ1+1)*γ1*ρ0^(γ1+1)+aL[5]*(γ2+1)*γ2*ρ0^(γ2+1)+aL[6]*(γ3+1)*γ3*ρ0^(γ3+1)+aL[7]*(γ4+1)*γ4*ρ0^(γ4+1) )
    K+=aL[8]*27*ρ0^2
    return K
end

function mstar_m(aL)
    ρ0=ρ0_GKW
    return 1/(1+2*mΛMeV/ħc^2*aL[2]*ρ0)
end

function OutputJLK()
    io1=open("JLK.csv","w")
    write(io1,"index,Parameter Name,J (MeV),L (MeV),K (MeV),m*/m,J-J_calc,L-L_calc,K-K_Calc\n")
    df=DataFrame(CSV.File("Lambda Parameters.csv"))

    for LParamType in 1:NumberOfParams
        #println(LParamType)
        pL=LambdaParameters.getParams(LParamType)
        aL=LambdaParameters.getaL(LParamType)

        J_anal=Calc_J(aL,pL.γ1,pL.γ2,pL.γ3,pL.γ4)
        L_anal=Calc_L(aL,pL.γ1,pL.γ2,pL.γ3,pL.γ4)
        K_anal=Calc_K(aL,pL.γ1,pL.γ2,pL.γ3,pL.γ4)

        ρ0=ρ0_GKW
        J_calc=Uopt(ρ0/2,ρ0/2,0,aL,pL.γ1,pL.γ2,pL.γ3,pL.γ4)

        h=0.001
        ρmesh=[ρ0-2*h,ρ0-h,ρ0,ρ0+h,ρ0+2*h]
        Umesh=zeros(Float64,5)
        for i in 1:5
            Umesh[i]=Uopt(ρmesh[i]/2,ρmesh[i]/2,0,aL,pL.γ1,pL.γ2,pL.γ3,pL.γ4)
        end
        L_calc=3*ρ0*MyLib.diff1st5pt(h,Umesh)
        K_calc=9*ρ0^2*MyLib.diff2nd5pt(h,Umesh)

        mstar_m_anal=mstar_m(aL)

        write(io1,"$(LParamType)")
        write(io1,",$(df[LParamType,"Parameter Name"])")
        write(io1,",$(J_anal)")
        write(io1,",$(L_anal)")
        write(io1,",$(K_anal)")
        write(io1,",$(mstar_m_anal)")
        write(io1,",$(J_anal-J_calc)")
        write(io1,",$(L_anal-L_calc)")
        write(io1,",$(K_anal-K_calc)\n")
    end

    close(io1)

end

OutputJLK()