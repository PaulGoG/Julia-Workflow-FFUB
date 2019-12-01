using DifferentialEquations
using Plots
using CSV
using Optim
using DiffEqParamEstim
using ODEInterfaceDiffEq

const G = 6.67408*1e-11        #Cosmological constant

function HdepZ(du, u, p, t)
        w = 0                      # u1,u2,u3  --> H,S,Mag
    du[1] = (1/((1+t)*u[1])) * (u[1]^2 + (8/3)*π*G*p[1]*(1 + 3*w)*(t+1)^(1 + 3*w) - p[2])          # p1,p2 -> 𝜺0, Λ
    du[2] = 1/((1+t)*u[1])
    u[3] = (-2.5)*p[3]*log10(1/(u[2]*u[2])) + p[4]       #p4, p5 --> m0,m1
end

Data = CSV.read("C:\\Users\\GoG\\Desktop\\JuliaCodes\\Dubna2019Cosmology\\GoG\\Data.csv"; copycols = true)
ok=1
Redshift = zeros(length(Data[3]))
Magnitude = zeros(length(Data[3]))
for i=1:length(Data[3])
    if ((!ismissing(Data[1][i])) || (!ismissing(Data[2][i])) || (!ismissing(Data[3][i])))
            try
                global Redshift[ok] = parse(Float64, Data[3][i])
                global Magnitude[ok] = parse(Float64, Data[2][i])
                global ok+=1
            catch Error
                println( )
            end

    end
end
ok=1
goodRedshift = zeros(6605-8)
goodMagnitude = zeros(6605-8)
for i =1:length(Redshift)
    if Redshift[i] > 0.
        goodRedshift[ok] = Redshift[i]
        goodMagnitude[ok] = Magnitude[i]
        global ok+=1
    end
end
ok=1
matrix = copy([goodRedshift,goodMagnitude])
aux = unique(matrix[1])
aux2 = indexin(aux, matrix[1])
for i in 1:length(goodRedshift)
    for j in 1:length(aux2)
        if i == aux2[j]
            matrix[2][ok] = matrix[2][i]
            global ok+=1
            break
        end
    end
end
for i in reverse(ok:length(goodRedshift))
        deleteat!(matrix[2], i)
end
unique!(matrix[1])
p = sortperm(matrix[1])
matrix[2] = copy(matrix[2][p])
sort!(matrix[1])
TimeIntervals = copy(matrix[1])
Mdata = copy(matrix[2])

function prob_func(prob, i, repeat)
    u0 = zeros(3)
    for j = 1:3
        u0[j] = prob.u0[j]*randn()
    end
    return ODEProblem(prob.f, u0 ,prob.tspan,prob.p)
end

f = HdepZ
u0 = [73., 1e-1, 0.1] # H0, S0, Mag
tspan = (1e-4, 1.)
p0 = [10., 1., 1., 15.] # ε0, Λ, m0, m1
N = 100

prob = ODEProblem(f, u0, tspan, p0)
monte_prob = MonteCarloProblem(prob, prob_func = prob_func)
#sim = solve(monte_prob, alg_hints=[:stiff], trajectories = N, saveat = TimeIntervals, abstol=1e-6, reltol=1e-6, maxiters=Int(1e8))

function loss(sim)
    x = 0.
    for i in 1:N
        for j in 1:length(sim[i][3,:])
            x = x +  (Mdata[j] - sim[i][3,j])^2
        end
    end
    return x
end

obj = build_loss_objective(monte_prob, Kvaerno5(), loss, trajectories=N , saveat = TimeIntervals, abstol=1e-7, reltol=1e-7, maxiters=Int(1e8))
result = optimize(obj, [1.,1.,1.,1.])
Optim.minimizer(result)



p0 = copy(Optim.minimizer(result))
prob = ODEProblem(f, u0, tspan, p0)
sol = solve(prob, saveat = TimeIntervals)
Magnitudini = [p0[3]*(-2.5)*log10(1/sol[i][2]^2) + p0[4] for i in 1:length(sol[2,:])]
plot(TimeIntervals, Magnitudini )
scatter!(TimeIntervals, Mdata)
