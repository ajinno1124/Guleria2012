include("Main.jl")
using .Threads

function run(NParamType,LParamType,io1)
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

		if ZN[i][1]==57 || ZN[i][1]==20
			Check=OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType,α=0.1)
		else
			Check=OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType)
		end
		#println(Check)
		if Check==false
			write(io1,"$(AN.Z)")
			write(io1,",$(AN.N)")
			write(io1,",$(AN.Λ)")
			write(io1,",$(NParamType)")
			write(io1,",$(LParamType)\n")
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
	#NParamType="SK3"
	#LParamType=[9,26,-1]

	#NParamType="SLy4"
	#LParamType=1:46

	#a3tuned params 2022/11/01
	#NParamType="SLy4"
	#LParamType=47:50

	#NParamType="SLy4"
	#LParamType=21:25

	#run all
	NParamType="SLy4"
	LParamType=1:25
	LParamType=vcat(LParamType,47:50)

	io1=open("NotConverge.csv","w")
	write(io1,"Z,N,L,NParamType,LParamType\n")

	@threads for i=eachindex(LParamType)
		run(NParamType,LParamType[i],io1)
	end

	close(io1)
end

function run_NotConverge()
	df=DataFrame(CSV.File("NotConverge.csv"))

	ListNotConverge=ones(Float64,nrow(df))

	@threads for i in 1:nrow(df)
		Z=df[i,"Z"]
		N=df[i,"N"]
		L=df[i,"L"]
		#NParamType=df[i,"NParamType"]
		NParamType="SLy4"
		LParamType=df[i,"LParamType"]
		AN=AtomNum(Z,N,L)

		println("\nZ=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")

		Check=OutPutFiles(AN,NParamType=NParamType,LParamType=LParamType,α=0.1,MaxIter=100)

		if Check==false
			ListNotConverge[i]=0
		end

	end

	for i in 1:nrow(df)
		println("$i,$(ListNotConverge[i])")
	end

end

@time run_threads()
@time run_NotConverge()
