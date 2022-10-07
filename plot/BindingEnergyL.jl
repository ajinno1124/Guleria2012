using CSV
using DelimitedFiles

function LYBindingEnergy()
	NParamType="SK3"
	LParamType="LY1"

	AN=[
		[6,5],
		[8,7],
		[14,13],
		[23,27],
		[39,49],
		[52,81],
		[82,125]
	]

	cd("../data/")
	io=open("BE_$(NparamType)$(LParamType).csv","w")
	write(io,"A, Z, N, lLam, jLam, B.E. Lambda(MeV)")
	for i=eachindex(AN)
		data1,header1=readdlm("Z$(AN[i][1])N$(AN[i][2])L0_$(NParamType)$(LParamType)/energy2.csv",',',header=true)
		df1=DataFrame(data1, vec(header1))
		data1,header1=readdlm("Z$(AN[i][1])N$(AN[i][2])L1_$(NParamType)$(LParamType)/energy.csv",',',header=true)
		df2=DataFrame(data2,vec(header2))
		for n in 1:length(df2["jlam"])
			write(io,"$(AN[i][1]+AN[i][2])")
			write(io,",$(AN[i][1])")
			write(io,",$(AN[i][2])")
			write(io,",$(df2["lLam"][n])")

end