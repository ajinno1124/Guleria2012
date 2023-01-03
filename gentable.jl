using Parameters
using CSV
using DataFrames
using Formatting

function WriteDoll(io1,val,digit)
    fmt="%.$(digit)f"
    s=sprintf1(fmt,val)
    write(io1,raw"$")
    if val==0
        write(io1,"0")
    else
        write(io1,s)
    end
    write(io1,raw"$")
end

function WriteAnd(io1)
    write(io1," & ")
end

function WriteBackSlash(io1)
    write(io1,"\\\\")
end

function WriteValue(io1,df,ParamNum::Vector{Int64},Key::String,digit)
    for i=eachindex(ParamNum)
        WriteAnd(io1)
        WriteDoll(io1,df[ParamNum[i],Key],digit)
    end
end

function OutputTable()
    io1=open("ParamTable.txt","w")
    write(io1,"\\hline")
    write(io1,"\n")
    write(io1,"\\hline")
    write(io1,"\n")
    write(io1,"& GKW2 & GKW3 & GKW2    & GKW3    & LY\${\\rm I}\$ & LY\${\\rm IV}\$ & HP\$\\Lambda\$2 \\\\")
    write(io1,"\n")
    write(io1,"&      &      & +Kohno2 & +Kohno3 &                &                 & \\\\")
    write(io1,"\n")
    write(io1,"\\hline")
    write(io1,"\n")
    df_param=DataFrame(CSV.File("Lambda Parameters.csv",comment="#"))
    ParamNum=[47,48,49,50,9,12,14]

    for i in 1:7
        if i==1
            write(io1,raw"$t_0^\Lambda~({\rm MeV~fm^3})$")
            key="u0"
            digit=1
        elseif i==2
            write(io1,raw"$t_1^\Lambda~({\rm MeV~fm^5})$")
            key="u1"
            digit=1
        elseif i==3
            write(io1,raw"$t_2^\Lambda~({\rm MeV~fm^5})$")
            key="u2"
            digit=1
        elseif i==4
            write(io1,raw"$t_{3,1}^\Lambda~({\rm MeV~fm}^4)$")
            key="u31"
            digit=1
        elseif i==5
            write(io1,raw"$t_{3,2}^\Lambda~({\rm MeV~fm}^5)$")
            key="u32"
            digit=1
        elseif i==6
            write(io1,raw"$x_0^\Lambda$")
            key="y0"
            digit=4
        elseif i==7
            write(io1,raw"$x_3^\Lambda$")
            key="y3"
            digit=4
        end
        WriteValue(io1,df_param,ParamNum,key,digit)
        WriteBackSlash(io1)
        write(io1,"\n")
    end


    df_chi=DataFrame(CSV.File("data/BindingEnergyLam/ChiSquared.csv",comment="#"))

    for i in 1:6
        if i==1
            write(io1,raw"$J_\Lambda$")
            key="J (MeV)"
            digit=2
        elseif i==2
            write(io1,raw"$L_\Lambda$")
            key="L (MeV)"
            digit=2
        elseif i==3
            write(io1,raw"$K_\Lambda$")
            key="K (MeV)"
            digit=1
        elseif i==4
            write(io1,raw"$m^*_\Lambda/m_\Lambda$")
            key="m*/m"
            digit=4
        elseif i==5
            write(io1,raw"MSD")
            key="ChiSquare6"
            digit=4
        elseif i==6
            write(io1,raw"\chi^2")
            key="ChiSquare1"
            digit=4
        end
        WriteValue(io1,df_chi,ParamNum,key,digit)
        WriteBackSlash(io1)
        write(io1,"\n")
    end

    write(io1,"\\hline")
    write(io1,"\n")
    write(io1,"\\hline")
    write(io1,"\n")


    close(io1)

end

@time OutputTable()