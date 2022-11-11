using Parameters
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

    ρ_cutoff=1.5
    k_cutoff=2.0
end


function getkfq(ρq)
    return (6*π^2*ρq/gq)^(1/3)
end

function getτq(ρq)
    kfq=getkfq(ρq)
    return 3.0/5.0*kfq^2*ρq
end

function getaL(J,L,K,ms_m)
	aL=zeros(Float64,5)
	aL[2]=1/ρ0*ħc^2/(2*mΛMeV)*(1/ms_m-1)
	τ0=getτq(ρ0/2)*2

	aL[1]=(10*(J-aL[2]*τ0) - 3*(L-5*aL[2]*τ0) + 0.5*(K-10*aL[2]*τ0))/ρ0
	aL[4]=(-15*(J-aL[2]*τ0) + 5*(L-5*aL[2]*τ0) - 1*(K-10*aL[2]*τ0))/ρ0^(4/3)
	aL[5]=(6*(J-aL[2]*τ0) - 2*(L-5*aL[2]*τ0) + 0.5*(K-10*aL[2]*τ0))/ρ0^(5/3)

	return aL

end

function GenerateJLK()
	io1=open("JLK_given.csv","w")
	J=[-33,-32,-31,-30,-29,-28,-27]
	L=[-50,-40,-30,-20,-10,0,10,20]
	K=[0,300,600]
	ms_m=[0.6,0.65,0.7,0.75,0.8,0.85,0.90,0.95,1.00]

	write(io1,"index,Parameter Name,paper")
	write(io1,",a1,a2,a3,a4,a5,a6,a7,a8")
	write(io1,",gamma1,gamma2,gamma3,gamma4")
	write(io1,",u0,u1,u2,u3,u31,u32")
	write(io1,",J (MeV),L (MeV),K (MeV),m*/m\n")

	count=1

	for j=eachindex(J)
		for l=eachindex(L)
			for k=eachindex(K)
				for mi=eachindex(ms_m)

					write(io1,"$(count),J$(J[j])L$(L[l])K$(K[k])ms_m$(ms_m[mi]),20221111")
					aL=getaL(J[j],L[l],K[k],ms_m[mi])
					write(io1,",$(aL[1]),$(aL[2]),$(aL[3]),$(aL[4]),$(aL[5])")
					write(io1,",0,0,0") #aL6,7,8
					write(io1,",$(1/3),$(2/3),0,0") #γ1,2,3,4
					write(io1,",$(aL[1]),$(aL[2]+2*aL[3]),$(3*aL[2]-2*aL[3]),0,$(8/3*aL[4]),$(8/3*aL[5])")
					write(io1,",$(J[j]),$(L[l]),$(K[k]),$(ms_m[mi])")
					write(io1,"\n")
					count+=1
				end
			end
		end
	end

	close(io1)

end

GenerateJLK()