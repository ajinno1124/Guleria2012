include("Main.jl")

function run_VB()
	NParamType="VB1"
	rm("data/$(NParamType)",force=true,recursive=true)
    mkpath("data/$(NParamType)")
    cd("data/$(NParamType)")
    ZN=[
        [8,8],
        [20,20],
		[20,28],
		[40,50],
		[82,126],
		[114,184]
    ]

    for i=eachindex(ZN)
        AN=AtomNum(ZN[i][1],ZN[i][2],0)
        println("\nZ=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")
		OutPutFiles(AN,NParamType=NParamType,LParamType="NaN",α=0.1)
    end

	cd("../../")

end

run_VB()