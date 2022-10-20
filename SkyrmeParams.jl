module NuclParameters
    mutable struct NuclParams
        t0::Float64
        t1::Float64
        t2::Float64
        t3::Float64
        x0::Float64
        x1::Float64
        x2::Float64
        x3::Float64
        σ::Float64
        W0::Float64
    end

    function getParams(ParamType::String)
        #SLy4
        if ParamType=="SLy4"
            t0=-2488.91
            t1=486.82
            t2=-546.39
            t3=13777.0
            x0=0.834
            x1=-0.344
            x2=-1.000
            x3=1.354
            σ=1/6
            W0=123.0
        elseif ParamType=="SLy5"
            t0=-2484.88
            t1=483.13
            t2=-549.40
            t3=13763.0
            x0=0.778
            x1=-0.328
            x2=-1.000
            x3=1.267
            σ=1/6
            W0=126.0
        elseif ParamType=="SLy6"
            t0=-2479.50
            t1=462.18
            t2=-448.61
            t3=13763.0
            x0=0.825
            x1=-0.465
            x2=-1.000
            x3=1.355
            σ=1/6
            W0=122.0
        elseif ParamType=="SLy7"
            t0=-2482.41
            t1=457.97
            t2=-419.85
            t3=13677.0
            x0=0.846
            x1=-0.511
            x2=-1.000
            x3=1.391
            σ=1/6
            W0=126.0
        elseif ParamType=="SkM*"
            t0=-2645.00
            t1=410.00
            t2=-135.00
            t3=15595.0
            x0=0.09
            x1=0.00
            x2=0.00
            x3=0.00
            σ=1/6
            W0=130.0
		elseif ParamType=="SK3"
            t0=-1128.75
            t1=395.0
            t2=-95.0
            t3=14000.0
            x0=0.45
            x1=0.00
            x2=0.00
            x3=1.0
            σ=1.0
            W0=120.0
		# exchange coulomb energyを含まないので、別で取り扱う
        #elseif ParamType=="VB1"
        #    t0=-1057.3
        #    t1=235.9
        #    t2=-100.0
        #    t3=14463.5
        #    x0=0.56
        #    x1=0.00
        #    x2=0.00
        #    x3=0.00
        #    σ=1
        #    W0=120.0
        elseif ParamType=="SKS3"
            t0=-2014.7
            t1=361.0
            t2=-29.5
            t3=12756
            x0=-0.319
            x1=0.732
            x2=4.95
            x3=-0.904
            σ=0.2604
            W0=94
        end

        return NuclParams(t0,t1,t2,t3,x0,x1,x2,x3,σ,W0)
    end

    function getaN(ParamType::String)
        aN=zeros(Float64,10)
        p=getParams(ParamType)
        aN[1]=0.25*p.t0*(2+p.x0)
        aN[2]=-0.25*p.t0*(2*p.x0+1)
        aN[3]=1/24*p.t3*(2+p.x3)
        aN[4]=-1/24*p.t3*(2*p.x3+1)
        aN[5]=1/8*(p.t1*(2+p.x1)+p.t2*(2+p.x2))
        aN[6]=1/8*(p.t2*(2*p.x2+1)-p.t1*(2*p.x1+1))
        aN[7]=1/32*(3*p.t1*(2+p.x1)-p.t2*(2+p.x2))
        aN[8]=-1/32*(3*p.t1*(2*p.x1+1)+p.t2*(2*p.x2+1))
        aN[9]=-1/16*(p.t1*p.x1+p.t2*p.x2)
        aN[10]=1/16*(p.t1-p.t2)

        #if ParamType=="VB1"
        #    aN[3]=p.t3/8
        #    aN[4]=-aN[3]
        #end

        return aN
    end
end


module LambdaParameters
    using CSV
    using DataFrames

    mutable struct LambdaParams
        γ1::Float64
		γ2::Float64
        γ3::Float64
        γ4::Float64
        u0::Float64
        u1::Float64
        u2::Float64
        u3::Float64
        u31::Float64
		u32::Float64
        u33::Float64
        u34::Float64
        y0::Float64
        y31::Float64
		y32::Float64
        y33::Float64
        y34::Float64
    end

    function getParams(ParamType::Int)
        #=
		γ1,γ2,u0,u1,u2,u3,u31,u32,y0,y31,y32=[0,0,0,0,0,0,0,0,0,0,0]

        if ParamType=="HPL1"
            γ1 = 1
            u0 = -326.395
            u1 = 72.627
            u2 = -8.584
            u3 = 0.0
            u31 = 1746.041
            y0 = -0.223
            y31 = -0.389
        elseif ParamType=="HPL2"
            γ1 = 1
            u0 = -399.946
            u1 = 83.426
            u2 = 11.455
            u3 = 0.0
            u31 = 2046.818
            y0 = -0.486
            y31 = -0.660
        elseif ParamType=="HPL3"
            γ1 =  1/3
            u0 = -498.515
            u1 = 65.203
            u2 = 19.001
            u3 = 0.0
            u31 = 995.832
            y0 = -0.492
            y31 = -0.444
        elseif ParamType=="HPL4"
            γ1 = 1/3
            u0 = -475.584
            u1 = 99.058
            u2 = -20.890
            u3 = 0.0
            u31 = 1375.172
            y0 = -0.350
            y31 = -0.724
        elseif ParamType=="NL1"
            γ1 = 1
            u0 = -253.3250
            u1 = 147.1264
            u2 = -83.5843
            u3 = 0.0
            u31 = 1684.9876
            y0 = 0.5802
            y31 = 0.4831
        elseif ParamType=="OL1"
            γ1 = 1
            u0 = -236.5835
            u1 = 116.8704
            u2 = -112.8812
            u3 = 0.0
            u31 = 1453.3493
            y0 = 0.1271
            y31 = -0.3110
        elseif ParamType=="NL2"
            γ1 = 1/3
            u0 = -518.620
            u1 = 82.0944
            u2 = -19.9772
            u3 = 0.0
            u31 = 1190.1894
            y0 = -0.1392
            y31 = 0.3126
        elseif ParamType=="OL2"
            γ1 = 1/3
            u0 = -417.7593
            u1 = 1.5460
            u2 = -3.2617
            u3 = 0.0
            u31 = 1102.2221
            y0 = -0.3854
            y31 = -0.5645
        elseif ParamType=="SKSH1"
            γ1 = 1.0
            u0 = -176.5
            u1 = -35.8
            u2 = 44.1
            u3 = 0.0
            u31 = 0.0
            y0 = 0.0
            y31 = 0.0
        elseif ParamType=="SKSH2"
            γ1 = 1.0
            u0 = -290.0
            u1 = 21.7
            u2 = -20.3
            u3 = 1850
            u31 = 0.0
            y0 = 0.0
            y31 = 0.0
		elseif ParamType=="LY1"
			γ1 = 1.0/3.0
            u0 = -476.0
            u1 = 42.0
            u2 = 23.0
            u3 = 0.0
            u31 = 1514.1
            y0 = -0.0452
            y31 = -0.280
        elseif ParamType=="GKW2"
            γ1 = 1/3
			γ2 = 2/3
            u0 = -968.125
            u1 = 0.0
            u2 = 0.0
            u3 = 0.0
            u31 = 4371.717378
			u32 = -1210.177854
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
        elseif ParamType=="GKW3"
            γ1 = 1/3
			γ2 = 2/3
            u0 = -500.625
            u1 = 0.0
            u2 = 0.0
            u3 = 0.0
            u31 = 4.912041998
			u32 = 2850.138497
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
		elseif ParamType=="GKW2+MD1"
			γ1 = 1/3
			γ2 = 2/3
            u0 = -0.904726223
            u1 = 32.29388211
            u2 = 96.88164632
            u3 = 0.0
            u31 = -3617.064524
			u32 = 4567.805539
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
		elseif ParamType=="GKW3+MD2"
			γ1 = 1/3
			γ2 = 2/3
            u0 = -90.81915637
            u1 = 35.98159235
            u2 = 107.944777
            u3 = 0.0
            u31 = -3455.520625
			u32 = 5194.926019
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
		elseif ParamType=="GKW3+MD3"
			γ1 = 1/3
			γ2 = 2/3
            u0 = -1.051619608
            u1 = 40.95591467
            u2 = 122.867744
            u3 = 0.0
            u31 = -4124.597901
			u32 = 5603.400881
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
		elseif ParamType=="GKW2_1.5"
			γ1 = 1/3
			γ2 = 2/3
            u0 = -352.2083744
            u1 = 0.0
            u2 = 0.0
            u3 = 0.0
            u31 = -951.8629975
			u32 = 3048.439943
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
		elseif ParamType=="GKW3_1.5"
			γ1 = 1/3
			γ2 = 2/3
            u0 = -388.2767987
            u1 = 0.0
            u2 = 0.0
            u3 = 0.0
            u31 = -1081.792893
			u32 = 3807.417485
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
		elseif ParamType=="GKW2_1.5+Kohno2"
			γ1 = 1/3
			γ2 = 2/3
            u0 = -352.2083744
            u1 = 40.84151607
            u2 = 122.5245482
            u3 = 0.0
            u31 = -951.8629973
			u32 = 2654.452384
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
		elseif ParamType=="GKW3_1.5+Kohno3"
			γ1 = 1/3
			γ2 = 2/3
            u0 = -388.2767987
            u1 = 53.53651325
            u2 = 160.6095397
            u3 = 0.0
            u31 = -1081.792893
			u32 = 3290.964568
            y0 = 0.0
            y31 = 0.0
			y32 = 0.0
		elseif ParamType=="LY1_a3zero"
			γ1 = 1.0/3.0
            u0 = -476.0
            u1 = 16.25
            u2 = 48.75
            u3 = 0.0
            u31 = 1514.1
            y0 = -0.0452
            y31 = -0.280
		elseif ParamType=="NaN"
			# all parameters are 0
		else
			println("$(ParamType) is not defined!")
        end
    =#

        if ParamType==-1
            return LambdaParams(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        else
            df=DataFrame(CSV.File("Lambda Parameters.csv"))
            γ1=df[ParamType,"Gamma1"]
            γ2=df[ParamType,"Gamma2"]
            γ3=df[ParamType,"Gamma3"]
            γ4=df[ParamType,"Gamma4"]
            u0=df[ParamType,"u0"]
            u1=df[ParamType,"u1"]
            u2=df[ParamType,"u2"]
            u3=df[ParamType,"u3"]
            u31=df[ParamType,"u31"]
            u32=df[ParamType,"u32"]
            u33=df[ParamType,"u33"]
            u34=df[ParamType,"u34"]
            y0=df[ParamType,"y0"]
            y31=df[ParamType,"y31"]
            y32=df[ParamType,"y32"]
            y33=df[ParamType,"y33"]
            y34=df[ParamType,"y34"]
            #println(LambdaParams(γ1,γ2,γ3,γ4,u0,u1,u2,u3,u31,u32,u33,u34,y0,y31,y32,y33,y34))
            return LambdaParams(γ1,γ2,γ3,γ4,u0,u1,u2,u3,u31,u32,u33,u34,y0,y31,y32,y33,y34)
        end
    end

    function getaL(ParamType::Int)
        aL=zeros(Float64,8)
        p=getParams(ParamType)
        aL[1]=p.u0*(1+0.5*p.y0)
        aL[2]=0.25*(p.u1+p.u2)
        aL[3]=1.0/8.0*(3*p.u1-p.u2)
        #if ParamType=="SKSH1" || ParamType=="SKSH2"
        #    aL[3]=0.25*(3*p.u1-p.u2)
        #end
        aL[4]=3.0/8.0*p.u31*(1+0.5*p.y31)
		aL[5]=3.0/8.0*p.u32*(1+0.5*p.y32)
        aL[6]=3.0/8.0*p.u33*(1+0.5*p.y33)
        aL[7]=3.0/8.0*p.u34*(1+0.5*p.y34)
        aL[8]=0.25*p.u3

        return aL
    end

end
