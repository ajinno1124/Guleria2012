include("Main.jl")

function run_LY()
    ZN=[
        [6,5],
        [8,7],
        [14,13],
        [23,27],
        [39,49],
        [52,81],
        [82,125]
    ]

    for i=eachindex(ZN)
        AN=AtomNum(ZN[i][1],ZN[i][2],0)
        println("Z=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")
        OutPutFiles(AN,NParamType="SK3",LParamType="NaN")

        AN=AtomNum(ZN[i][1],ZN[i][2],1)
        println("Z=$(AN.Z), N=$(AN.N), L=$(AN.Λ)")
        OutPutFiles(AN,NParamType="SK3",LParamType="LY1")
    end
end

run_LY()