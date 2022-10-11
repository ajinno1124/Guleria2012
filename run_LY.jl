include("Main.jl")

function run_LY()
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
			OutPutFiles(AN,NParamType="SK3",LParamType="NaN",α=0.1)
		else
        	OutPutFiles(AN,NParamType="SK3",LParamType="NaN")
		end


        AN=AtomNum(ZN[i][1],ZN[i][2],1)
        println("Z=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")
		if ZN[i][1]==57
			OutPutFiles(AN,NParamType="SK3",LParamType="LY1",α=0.1)
		else
        	OutPutFiles(AN,NParamType="SK3",LParamType="LY1")
		end
    end

end

run_LY()