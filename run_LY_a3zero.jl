include("Main.jl")

function run_LY()
	NParamType="SK3"
	LParamType=26
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

    for i=eachindex(ZN)
        AN=AtomNum(ZN[i][1],ZN[i][2],0)
        println("\nZ=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")
		if ZN[i][1]==57
			OutPutFiles(AN,NParamType=NParamType,LParamType=-1,α=0.1)
		else
        	OutPutFiles(AN,NParamType=NParamType,LParamType=-1)
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

run_LY()