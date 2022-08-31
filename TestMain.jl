include("Main.jl")
using Plots

AN=AtomNum(82,126,1) #lead
#AN=AtomNum(8,8,1) #oxigen

function TestInitPot()
    rmesh=getrmesh()
    h2m,dh2m,V,W=InitPot(AN,rmesh)
    plot(xlabel="r",ylabel="V",title="Woods Saxon Z=$(AN.Z), N=$(AN.N)")
    plot!(rmesh,V[1,:],label="proton")
    plot!(rmesh,V[2,:],label="neutron")
    plot!(rmesh,V[3,:],label="Λ")
    plot!(rmesh,(@. W[1,:]/rmesh[:]),label="W/r")
end

function CheckABC()
    rmesh=getrmesh()
    h2m,dh2m,V,W=InitPot(AN,rmesh)

    h2mB=h2m[1,:]
    dh2mB=dh2m[1,:]
    VB=V[1,:]
    WB=W[1,:]

    QN=QuantumNumber(0.5,0,1) #(j,l,B)
    j=QN.j
    l=QN.l

    A=zeros(Float64,Nmesh)
    B=zeros(Float64,Nmesh)
    C=zeros(Float64,Nmesh)
    @. A[:] += -h2mB[:]
    @. B[:] += -dh2mB[:]
    @. C[:] += h2mB[:]*l*(l+1)/(rmesh[:]^2) + VB[:] + dh2mB[:]/rmesh[:] + WB[:]/rmesh[:]*(j*(j+1)-l*(l+1)-0.75)
    
    plot(xlabel="r",ylabel="A,B,C",ylim=(minimum(C)-5,1),title="(j,l,B)=($j,$l,$(QN.B)), Z=$(AN.Z), N=$(AN.N)")
    plot!(rmesh,A,label="A")
    plot!(rmesh,B,label="B")
    plot!(rmesh,C,label="C")
end

function Testgetrmesh(rc)
    h=rc/(Nmesh-0.5)
    rmesh=range(0.5*h,(Nmesh-0.5)*h,length=Nmesh)
    return rmesh
end

function CheckHmat()
    rc=[30,60,120]
    plot(xlabel="Index of Energy", ylabel="Eigen Energy")
    for i in eachindex(rc)
        rmesh=Testgetrmesh(rc[i])
        #println(typeof(rmesh))

        A=zeros(Float64,Nmesh)
        B=zeros(Float64,Nmesh)
        C=zeros(Float64,Nmesh)
        @. A[:] += -0.5
        @. B[:] += -1/rmesh[:]
        @. C[:] += -1/rmesh[:]
        Hmat=MakeHmat(A,B,C,rmesh)
        #=
            Hmat=Hmat+zeros(Float64,Nmesh,Nmesh)
            E,ψ=eigen(Hmat)
        =#
        
        for i in 1:Nmesh
            Hmat[i,i]+=offset
        end
        E,ψ=eigs(Hmat,nev=10,maxiter=1000,which=:SR)
        E=real(E.-offset)
        
        #E=real(sort(real(E)))
        plot!(real(E[1:10]),label="rc=$(rc[i])")
    end
    plot!()
end


function TestInitialCondition()
    
    InitState=InitialCondition(AN)
    for i=eachindex(InitState[1])
        l=InitState[1][i].QN.l
        j=InitState[1][i].QN.j
        E=InitState[1][i].E
        println("(l,j,E)=($l,$j,$E)")
    end
    
    #proton波動関数をプロット
    plot(xlabel="r", ylabel="R", title="Initial Proton Wave Function")
    rmesh=getrmesh()
    for i=eachindex(InitState[1])
    #for i=1:2
        if InitState[1][i].E<0
            l=InitState[1][i].QN.l
            j=InitState[1][i].QN.j
            plot!(rmesh,InitState[1][i].ψ,label="(l,j)=($l,$j)")
        end
    end
    #println(InitState[1][1].ψ)
    plot!()
end