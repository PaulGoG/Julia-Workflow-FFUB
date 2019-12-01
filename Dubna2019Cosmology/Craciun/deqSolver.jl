using DifferentialEquations

#defining the set of differetnial equations
function deq(du,u,p,t)
    W = [-1,-1/2,0,1/3]
    H = u[1]
    z = t
    e0 = p

    if (e0[1]+e0[2])>2*u0[1]
        e0[1]=u0[1]
        e0[2]=u0[1]
    end

    #-1 and 0
    #du[1] = H/(z+1)+1/H*(e0[1]*(1.0+3.0*W[1])*(z+1.0)^(2.0+3.0*W[1])+e0[2]*(1.0+3.0*W[3])*(z+1.0)^(2.0+3.0*W[3]))
    #-1/2 and 1/3
    du[1] = H/(z+1)+1/H*(e0[1]*(1.0+3.0*W[2])*(z+1.0)^(2.0+3.0*W[2])+e0[2]*(1.0+3.0*W[4])*(z+1.0)^(2.0+3.0*W[4]))

    du[2] = 1/(H*(z+1))
end

#calculating magnitude as a function of distance
function mag(S,m0,m1)
    return -2.5.*m0.*log10.(1.0./S.^2).+m1
end

#import the data
using CSV
#set good path
cd("C:\\Users\\GoG\\Desktop\\JuliaCodes\\Dubna2019Cosmology\\Craciun")
data = CSV.read("data\\prelucratedData.csv")

#cost function
function cost(param)
    #create the ODEProblem
    prob = ODEProblem(deq,u0,tspan,[param[1],param[2],param[3],param[4]])
    #solve the ODEProblem
    sol = solve(prob,Rodas5())

    try
        M = mag(sol.(data[2],idxs=2),param[5],param[6])
        return sum((M.-data[1]).^2)
    catch
        return 10e200
    end

end

#initial values
u0 = [73.0,0.0]
#span to evaluate
tspan = (0.0,40.0)
#initial values for parameters
initialp = [10,10,1.0,36.0]


#here we minimize the cost function to find the parameters
using Optim
#setting bounds
newparam = optimize(cost,initialp,iterations=100,time_limit=3)
print(newparam)

#we construct a new solution with the new parameters
prob1 = ODEProblem(deq,u0,tspan,[newparam.minimizer[1],newparam.minimizer[2]])#newparam.minimizer[1])
sol1 = solve(prob1,Rodas5())
Magnitude = mag(sol1.(data[2],idxs=2),newparam.minimizer[3],newparam.minimizer[4])

#plot the new solution against the data
using Plots
plot(sol1.t,sol1[1,:])

#=
scatter(data[2],data[1])
plot!(data[2],Magnitude)
plot(data[2],sol1.(data[2],idxs=1))
=#
