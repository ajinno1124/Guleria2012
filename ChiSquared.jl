using CSV, DataFrames

function Make_BindingEnergyLamData(NParamType,LParamType)
    ZN=[
		#[2,5],
		#[3,5],
		#[4,5],
		#[5,4],
		#[5,5],
		#[5,6],
        #[6,5],

		[6,6],
		[6,7],
		[7,8],# not included in Ohnishi-san's File
        [8,7],
        [14,13],
		[16,15],
		[20,19],
        [23,27],
		[26,29],
        [39,49],
        [57,81],
        [82,125]
	]

    io1=open("data/BindingEnergyLam/BindingEnergy$(LParamType).csv","w")
    write(io1,"#Elcheck = e_Lam - (e_Lam using Rearrangement Energy)\n")
    write(io1,"A,Z,N,jLam,lLam,B.E Lambda(MeV),EL_check(MeV)\n")

    for i=eachindex(ZN)
        df1=DataFrame(CSV.File("data/$(NParamType)NaN/Z$(ZN[i][1])N$(ZN[i][2])L0_$(NParamType)NaN/Energy2.csv",comment="#"))
        df2=DataFrame(CSV.File("data/$(NParamType)$(LParamType)/Z$(ZN[i][1])N$(ZN[i][2])L1_$(NParamType)$(LParamType)/Energy.csv",comment="#"))
        l=0
        for n in 1:nrow(df2)
            if df2[n,"lLam"]==l
                write(io1,"$(ZN[i][1]+ZN[i][2])")
                write(io1,",$(ZN[i][1])")
                write(io1,",$(ZN[i][2])")
                write(io1,",$(df2[n,"jLam"])")
                write(io1,",$(df2[n,"lLam"])")
                write(io1,",$(-df2[n,"Etot(MeV)"]+df1[1,"Etot2(MeV)"])")
                write(io1,",$(df2[n,"El_Check(MeV)"])\n")
                l+=1
            end
        end
    end

    close(io1)

end

#1/Nd*\sum(B_exp-B_th)^2/\sigma^2
function ChiSquared1(B_exp::Vector{Float64},B_exp_error::Vector{Float64},B_th::Vector{Float64})
	χ1=0.0
	count=0.0
	for i in 1:length(B_exp)
		if isnan(B_th[i])==false
			χ1+=(B_exp[i]-B_th[i])^2/B_exp_error[i]^2
			count+=1.0
		end
	end
	χ1/=count
	return χ1
end

#1/Nd*\sum(B_exp-B_th)^2)/1MeV^2
function ChiSquared2(B_exp::Vector{Float64},B_th::Vector{Float64})
	χ2=0.0
	count=0.0
	for i in 1:length(B_exp)
		if isnan(B_th[i])==false
			χ2+=(B_exp[i]-B_th[i])^2
			count+=1.0
		end
	end
	χ2/=count
	χ2=χ2^0.5
	return χ2
end
#\sum(B_exp-B_th)^2/B_exp^2
function ChiSquared3(B_exp::Vector{Float64},B_th::Vector{Float64})
	χ3=0.0
	count=0.0
	for i in 1:length(B_exp)
		if isnan(B_th[i])==false
			χ3+=(B_exp[i]-B_th[i])^2/B_exp[i]^2
			count+=1.0
		end
	end
	χ3/=count
	return χ3,count
end

function FindBELam(Z,N,lLam,df_th)
	B_th=0.0
	for n in 1:nrow(df_th)
		if [Z,N,lLam]==[df_th[n,"Z"],df_th[n,"N"],df_th[n,"lLam"]]
			#if Z!=57
				B_th=df_th[n,"B.E Lambda(MeV)"]
				break
			#end
		end
	end
	if B_th==0.0
		B_th=NaN
		println("No Binding Energy of Lambda was found.")
	end

	return B_th
end

function LineupBindingEnergy(df_exp,df_th)
	B_exp=df_exp[:,"B. E. (MeV)"]
	B_exp_error=df_exp[:,"error(MeV)"]
	B_th=zeros(Float64,length(B_exp))

	for i in 1:nrow(df_exp)
		Z=df_exp[i,"Core Z"]
		N=df_exp[i,"Core N"]
		lLam=df_exp[i,"lLam"]

		B_th[i]=FindBELam(Z,N,lLam,df_th)

	end

	return B_exp,B_exp_error,B_th
end

function BE_swave(df_exp,df_th)
	B_exp_s=Float64[]
	B_exp_error_s=Float64[]
	B_th_s=Float64[]
	for i in 1:nrow(df_exp)
		if df_exp[i,"lLam"]==0
			push!(B_exp_s,df_exp[i,"B. E. (MeV)"])
			push!(B_exp_error_s,df_exp[i,"error(MeV)"])

			Z=df_exp[i,"Core Z"]
			N=df_exp[i,"Core N"]
			lLam=df_exp[i,"lLam"]

			push!(B_th_s,FindBELam(Z,N,lLam,df_th))
		end
	end

	return B_exp_s,B_exp_error_s,B_th_s

end

function BE_sp_splitting(df_exp,df_th)
	B_exp_sp=Float64[]
	B_exp_error_sp=Float64[]
	B_th_sp=Float64[]
	for i in 1:nrow(df_exp)
		if df_exp[i,"lLam"]==0
			Z=df_exp[i,"Core Z"]
			N=df_exp[i,"Core N"]

			for k in 1:nrow(df_exp)
				if [Z,N,1]==[df_exp[k,"Core Z"],df_exp[k,"Core N"],df_exp[k,"lLam"]]
					push!(B_exp_sp,df_exp[i,"B. E. (MeV)"]-df_exp[k,"B. E. (MeV)"])
					push!(B_exp_error_sp,(df_exp[i,"error(MeV)"]^2+df_exp[k,"error(MeV)"]^2)^0.5)

					push!(B_th_sp,FindBELam(Z,N,0,df_th)-FindBELam(Z,N,1,df_th))
				end
			end
		end
	end

	return B_exp_sp,B_exp_error_sp,B_th_sp
end

function BE_Heavier13C(df_exp,df_th)
	B_exp_13C=Float64[]
	B_exp_error_13C=Float64[]
	B_th_13C=Float64[]
	for i in 1:nrow(df_exp)
		if df_exp[i,"Core A"]>13
			push!(B_exp_13C,df_exp[i,"B. E. (MeV)"])
			push!(B_exp_error_13C,df_exp[i,"error(MeV)"])

			Z=df_exp[i,"Core Z"]
			N=df_exp[i,"Core N"]
			lLam=df_exp[i,"lLam"]

			push!(B_th_13C,FindBELam(Z,N,lLam,df_th))
		end
	end

	return B_exp_13C,B_exp_error_13C,B_th_13C
end

function Calc_ChiSquared(df_exp,df_th)
	B_exp,B_exp_error,B_th=LineupBindingEnergy(df_exp,df_th)
	B_exp_s,B_exp_error_s,B_th_s=BE_swave(df_exp,df_th)
	B_exp_sp,B_exp_error_sp,B_th_sp=BE_sp_splitting(df_exp,df_th)
	B_exp_13C,B_exp_error_13C,B_th_13C=BE_Heavier13C(df_exp,df_th)
	χ1=ChiSquared1(B_exp,B_exp_error,B_th)
	χ2=ChiSquared2(B_exp,B_th)
	χ3,count=ChiSquared3(B_exp,B_th)
	χ4=ChiSquared2(B_exp_s,B_th_s)
	χ5=ChiSquared2(B_exp_sp,B_th_sp)
	χ6=ChiSquared2(B_exp_13C,B_th_13C)

	return χ1,χ2,χ3,χ4,χ5,χ6,count
end

function Make_ChiSquaredData(index,df_ParamType,df_JLK)
    io1=open("data/BindingEnergyLam/ChiSquared.csv","w")

	write(io1,"#ChiSquare1 = 1/Nd*sum(B_exp-B_th)^2/sigma\n")
	write(io1,"#ChiSquare2 = sqrt(1/Nd*sum(B_exp-B_th)^2)^2) MeV\n")
	write(io1,"#ChiSquare3 = sum(B_exp-B_th)^2/B_exp^2\n")
	write(io1,"#ChiSquare4 = ChiSquare2 using only s-wave")
	write(io1,"#ChiSquare5 = ChiSquare2 using only s-p splitting\n")
	write(io1,"#ChiSquare6 = ChiSquare2 using only heavyer than 13C_Lam\n")
	write(io1,"index,NParamType,LParameterType,J (MeV),L (MeV),K (MeV),m*/m,ChiSquare1,ChiSquare2,ChiSquare3,ChiSquare4,ChiSquare5,ChiSquare6,Number of Data\n")

	df_exp=DataFrame(CSV.File("LamBindingEnergy.csv",comment="#"))
    for i=eachindex(index)
		LParamType_str=df_ParamType[index[i],"Parameter Name"]
		df_th=DataFrame(CSV.File("data/BindingEnergyLam/BindingEnergy$(LParamType_str).csv",comment="#"))
		χ1,χ2,χ3,χ4,χ5,χ6,count=Calc_ChiSquared(df_exp,df_th)
		write(io1,"$(index[i])")
		write(io1,",SLy4") #20221105 only use SLy4
		write(io1,",$(df_ParamType[index[i],"Parameter Name"])")
		write(io1,",$(df_JLK[index[i],"J (MeV)"])")
		write(io1,",$(df_JLK[index[i],"L (MeV)"])")
		write(io1,",$(df_JLK[index[i],"K (MeV)"])")
		write(io1,",$(df_JLK[index[i],"m*/m"])")
		write(io1,",$(χ1)")
		write(io1,",$(χ2)")
		write(io1,",$(χ3)")
		write(io1,",$(χ4)")
		write(io1,",$(χ5)")
		write(io1,",$(χ6)")
		write(io1,",$(count)\n")
    end

	close(io1)
end

function ExecuteAll()
    rm("data/BindingEnergyLam",force=true,recursive=true)
    mkpath("data/BindingEnergyLam")

    df=DataFrame(CSV.File("Lambda Parameters.csv"))
    #index=1:50
	#index=vcat(1:25,47:50)
	#index=vcat(1:20,47:50)
	#index=1:1562
	index=1:3578
    for i in index
        LParamType_str=df[i,"Parameter Name"]
        Make_BindingEnergyLamData("SLy4",LParamType_str)
    end
	df2=DataFrame(CSV.File("JLK.csv"))
    Make_ChiSquaredData(index,df,df2)
end

@time ExecuteAll()