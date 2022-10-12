include("Main.jl")

function run_SLy4()
	NParamType="SLy4"
	LParamtype="NaN"
	rm("data/$(NParamType)$(LParamtype)",force=true,recursive=true)
    mkpath("data/$(NParamType)$(LParamtype)")
    cd("data/$(NParamType)$(LParamtype)")
    ZN=[
        [8,8],
        [20,20],
		[20,28],
		[28,28],
		[28,50],
		[50,50],
		[50,82],
		[82,126],
    ]

    for i=eachindex(ZN)
        AN=AtomNum(ZN[i][1],ZN[i][2],0)
        println("\nZ=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")
		OutPutFiles(AN,NParamType=NParamType,LParamType="NaN")
    end

	cd("../..")

end

run_SLy4()