module MyLib
    #integral using trapezoid formula
    function IntTrap(x,y)
        N=length(x)
        ans=0.0
        ans+=2*sum(y)-y[N]-y[1]
        ans*=(x[N]-x[1])/(2*(N-1))
        return ans
    end

    #return the Vector of dydx
    function diff1st(h::Float64,y::AbstractArray)
        N=length(y)
        dydx=zeros(Float64,N)
        #dydx[1]=(y[2]-y[1])/h # order(h)
        dydx[1]=(-y[3]+4*y[2]-3*y[1])/(2*h) #order(h^2)
        for i in 2:N-1
            dydx[i]=(y[i+1]-y[i-1])/(2*h)
        end
        #dydx[N]=(y[N]-y[N-1])/h
        dydx[N]=(-y[N-2]+4*y[N-1]-3*y[N])/(2*(-h))

        return dydx
    end

    function diff2nd(h::Float64,y::AbstractArray)
        N=length(y)
        ddyddx=zeros(Float64,N)
        ddyddx[1]=(2*y[1]-5*y[2]+4*y[3]-y[4])/(h^2)
        for i in 2:N-1
            ddyddx[i]=(y[i+1]-2*y[i]+y[i-1])/(h^2)
        end
        ddyddx[N]=(2*y[N]-5*y[N-1]+4*y[N-2]-y[N-3])/(h^2)

        return ddyddx
    end

end

# ref. 計算物理学
# solve Poisson eq. U''=-4πrρ(r)
# qmax = \int d^3r ρ(r)
function SolvePoissonEq(ρ::Vector{Float64},rmesh,qmax)
    U=zeros(Float64,Nmesh)
    h=rmesh[2]-rmesh[1]
    U[1]=rmesh[1]
    U[2]=3*U[1]-h^2*4*π*rmesh[2]*ρ[2]
    #U[2]=U[1]-h^2*4*π*rmesh[2]*ρ[2]
    for i in 2:Nmesh-1
        U[i+1]=2*U[i]-U[i-1]-h^2*4*π*rmesh[i]*ρ[i]
    end

    #adjust the value by adding the solution of U''=0 (U=αr)
    α=(U[Nmesh]-qmax)/rmesh[Nmesh]
    @. U[:]-=α*rmesh[:]

    return U
end