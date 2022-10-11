# CSV, DelimitedFilesのdocumentationが読みずら過ぎるので、もう扱わん！
using CSV
using DataFrames
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

	io=open("../data/BE_$(NParamType)$(LParamType).csv","w")
	write(io,"A, Z, N, lLam, jLam, B.E. Lambda(MeV)")
	for i=eachindex(AN)
		#d1,h1=readdlm("../data/Z$(AN[i][1])N$(AN[i][2])L0_$(NParamType)NaN/energy2.csv",',',header=true)
		df1=CSV.read("../data/Z$(AN[i][1])N$(AN[i][2])L0_$(NParamType)NaN/energy2.csv",comment="#",DataFrame)
		#d2,h2=readdlm("../data/Z$(AN[i][1])N$(AN[i][2])L1_$(NParamType)$(LParamType)/energy.csv",',',header=true,makeunique=true)
		df2=CSV.read("../data/Z$(AN[i][1])N$(AN[i][2])L1_$(NParamType)$(LParamType)/energy.csv",comment="#",DataFrame)
		print(df2[!,"jLam"])
		for n in 1:length(df2["jLam",:])
			write(io,"$(AN[i][1]+AN[i][2]+1)") #single lambda hyper
			write(io,",$(AN[i][1])")
			write(io,",$(AN[i][2])")
			write(io,",$(df2["lLam",n])")
			write(io,",$(df2["jlam",n])")
			write(io,",$(df2["Etot",n]-df1["Etot",n])")
		end
	end

	close(io)

end