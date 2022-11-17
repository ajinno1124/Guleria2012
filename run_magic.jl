include("Main.jl")
using .Threads

function run(NParamType,LParamType)
    ZN=[
		[2,2],
		[8,8],
		[20,20],
		[20,28], #> 5.8E22 y
		[28,28], #6.075 d 10
		[28,50], #122.2 ms 51
		[50,50], #1.16 s 20
		[50,82], #39.7 s 8
		[82,126], #Stable

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

    for i=eachindex(ZN)
		if LParamType==-1
			AN=AtomNum(ZN[i][1],ZN[i][2],0)
		else
			AN=AtomNum(ZN[i][1],ZN[i][2],1)
		end

		println("\nZ=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")

		Check=OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType,α=0.1,MaxIter=50)
		@assert Check==true

		#=
		if Check==false
			write(io1,"$(AN.Z)")
			write(io1,",$(AN.N)")
			write(io1,",$(AN.Λ)")
			write(io1,",$(NParamType)")
			write(io1,",$(LParamType)\n")
		end
		=#
    end

end

function run_magic()
	#Lanskoy and Yamamoto
	NParamType="SLy4"
	LParamType=vcat(-1,9:12)

	@threads for i=eachindex(LParamType)
		run(NParamType,LParamType[i])
	end

	NParamType="SK3"
	LParamType=vcat(-1,9:12)

	@threads for i=eachindex(LParamType)
		run(NParamType,LParamType[i])
	end

end

#=
function run_NotConverge()
	df=DataFrame(CSV.File("NotConverge.csv"))

	ListNotConverge=zeros(Float64,nrow(df))

	@threads for i in 1:nrow(df)
		Z=df[i,"Z"]
		N=df[i,"N"]
		L=df[i,"L"]
		#NParamType=df[i,"NParamType"]
		NParamType="SLy4"
		LParamType=df[i,"LParamType"]
		AN=AtomNum(Z,N,L)

		println("\nZ=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")

		α=[0.09,0.08,0.07,0.06]
		for j=eachindex(α)
			Check=OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType,α=α[j],MaxIter=100)
			if Check==true
				ListNotConverge[i]=1.0
				break
			end
		end

	end

	for i in 1:nrow(df)
		println("$i,$(ListNotConverge[i])")
	end

end
=#

@time run_magic()
#@time run_NotConverge()
