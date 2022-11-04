using CSV, DataFrames

function Make_BindingEnergyLamData(NParamType,LParamType)
    ZN=[
		[2,5],
		[3,5],
		[4,5],
		[5,4],
		[5,5],
		[5,6],
        [6,5],
		[6,6],
		[7,8],
        [8,7],
        [14,13],
		[16,15],
		[20,19],
        [23,27],
        [39,49],
        [57,81],
        [82,125]
	]

    io1=open("data/BindingEnergyLam/BindingEnergy$(LParamType).csv","w")
    write(io1,"#Elcheck = e_Lam - (e_Lam using Rearrangement Energy)\n")
    write(io1,"A,Z,N,jLam,lLam,B.E Lambda(MeV),EL_check(MeV)\n")

    for i=eachindex(ZN)
        df1=DataFrame(CSV.File("data/Z$(ZN[i][1])N$(ZN[i][2])L0_$(NParamType)NaN/Energy2.csv",comment="#"))
        df2=DataFrame(CSV.File("data/Z$(ZN[i][1])N$(ZN[i][2])L1_$(NParamType)$(LParamType)/Energy.csv",comment="#"))
        l=0
        for n in 1:nrow(df2)
            if df2[n,"lLam"]==l
                write(io1,"$(ZN[i][1]+ZN[i][2])")
                write(io1,",$(ZN[i][1])")
                write(io1,",$(ZN[i][2])")
                write(io1,",$(df2[n,"jLam"])")
                write(io1,",$(df2[n,"lLam"])")
                write(io1,",$(df2[n,"jLam"])")
                write(io1,",$(-df2[n,"Etot(MeV)"]+df1[1,"Etot2(MeV)"])")
                write(io1,",$(df2[n,"El_Check(MeV)"])\n")
                l+=1
            end
        end
    end

    close(io1)
    
end

#1/Nd*\sum(B_exp-B_th)^2/\sigma
function ChiSquared1()
end

#1/Nd*\sum(B_exp-B_th)^2)/1MeV^2
function ChiSquared2()

end
#\sum(B_exp-B_th)^2/B_exp^2
function ChiSquared3()
end

function Make_ChiSquaredData(index::Array{1,Int64})
    io1=open("data/BindingEnergyLam/ChiSquared.csv","w")

    for i in index
        
    end


end

function ExecuteAll()
    rm("data/BindingEnergyLam",force=true,recursive=true)
    mkpath("data/BindingEnergyLam")

    df=DataFrame(CSV.File("Lambda Parameters.csv"))
    index=1:50
    for i in index
        LParamType_str=df[i,"Parameter Name"]
        Make_BindingEnergyLamData("SLy4",LParamType_str)
    end
    Make_ChiSquaredData(index)
end

ExecuteAll()