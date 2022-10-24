include("Main.jl")

function run(NParamType,LParamType)
    ZN=[
        [6,5],
        [8,7],
        [14,13],
        [23,27],
        [39,49],
        [57,81],
        [82,125]
    ]

    for i=eachindex(ZN)
        AN=AtomNum(ZN[i][1],ZN[i][2],0)
        println("\nZ=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")
		if ZN[i][1]==57
			OutPutFiles(AN,NParamType=NParamType,LParamType="NaN",α=0.05)
		else
        	OutPutFiles(AN,NParamType=NParamType,LParamType="NaN")
		end


        AN=AtomNum(ZN[i][1],ZN[i][2],1)
        println("Z=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")
		if ZN[i][1]==57
			OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType,α=0.1)
		else
        	OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType)
		end
    end

end

run("SK3","LY1")
run("SLy4","HPL2")
run("SLy4","GKW2_1.5")
run("SLy4","GKW3_1.5")
run("SLy4","GKW2_1.5+Kohno2")
run("SLy4","GKW3_1.5+Kohno3")