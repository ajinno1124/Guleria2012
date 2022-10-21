include("Main.jl")
using .Threads

function run(NParamType,LParamType)
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
		if LParamType==-1
			AN=AtomNum(ZN[i][1],ZN[i][2],0)
		else
			AN=AtomNum(ZN[i][1],ZN[i][2],1)
		end

		println("\nZ=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")

		if ZN[i][1]==57
			OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType,α=0.1)
		else
			OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType)
		end
    end

end

function run_threads()
	#GKW_w/o mom. dep.
	#NParamType="SLy4"
	#LParamType=[38,37,36,33,-1]

	#GKW w/ mom. dep.
	#NParamType="SLy4"
	#LParamType=[35,39,40,41,-1]

	#Lanskoy and Yamamoto
	NParamType="SK3"
	LParamType=[9,26,-1]
	@threads for i=eachindex(LParamType)
		run(NParamType,LParamType[i])
	end
end

run_threads()
#=
run("SLy4",38)
run("SLy4",37)
run("SLy4",36)
run("SLy4",33)
run("SLy4",-1)
=#