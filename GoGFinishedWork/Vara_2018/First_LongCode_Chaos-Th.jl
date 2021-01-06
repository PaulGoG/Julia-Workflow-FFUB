# This is very old, messy, undocumented code!!!

using Plots, DifferentialEquations, LsqFit, DynamicalSystems, QuadGK
"""
    Lorenz(du, u, p, t)
Creates the Lorenz attractor described by its system of differential equation

## Variables and parameters

- `du`: derivatives of the system eq of motion
- `u`: system of eq of motion
- `p`: constant parameters
- `t`: timestamps
"""
function Lorenz(du, u, p, t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3])-u[2]
    du[3] = u[1]*u[2]-p[3]*u[3]
end
"""
    ConjoinedEq(f, p, u0::Array, d0)

Creates a system of 2*n equations for evaluating 2 simultaneous trajectories at the same time steps

## Variables and parameters

- `f`: the dynamical system
- `p`: constant parameters
- `u0`: initial conditions
- `d0`: initial distance between the 2 trajectories
"""
function ConjoinedEq(f, p, u0::Array, d0)
    n = length(u0)
    function SystemEq(du ,u ,p ,t)
        f(view(du,1:n),view(u,1:n),p,t)
        f(view(du,n+1:2n),view(u,n+1:2n),p,t)
    end
    u1 = u0 + d0/√n
    return SystemEq, vcat(u0,u1)
end

function FunctionRepeater(f, p, u0, d0, tspan)
	SolveEq, u00 = ConjoinedEq(f,p,u0,d0)
	ODEProblem(SolveEq, u00, tspan, p)
end

function output_func(solution, i)
    d = DiffEqArray([norm(solution[1:3,i]-solution[4:6,i]) for i in 1:length(solution.t)],solution.t)
    return d, false
end

function fit_λ(sol, T_λ)
    dmed = timeseries_steps_mean(sol)
    model_func(x, λ) = d0 *exp.(λ[1]*x)
    c_fit = curve_fit(model_func, dmed.t[dmed.t .< T_λ], dmed[dmed.t .< T_λ], [0.])
    λ_fit = c_fit.param

    return λ_fit[1]
end

function initial_conditions(f, p, n, d0, t=1e3, ttr=100.)
    u0 = [0, 10., 0]
    tspan = (0., ttr)
    transientprob = ODEProblem(f, u0, tspan, p)
    sol = solve(transientprob, Vern9(), abstol=1e-9, reltol=1e-3, save_everystep=false)
    tspan = (ttr, ttr+t)
    u0 = sol[end]
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, Vern9(), abstol=1e-9, reltol=1e-3)

    return [sol(rand()*t+ttr) for i=1:n]
end

function calculate_miu(f, p, n, d0, t=1e3, ttr=100.)
    u0 = [0, 10., 0]
    tspan = (0., ttr)
    transientprob = ODEProblem(f, u0, tspan, p)
    sol = solve(transientprob, Vern9(), abstol=1e-9, reltol=1e-3, save_everystep=false)
    tspan = (ttr, ttr+t)
    u0 = sol[end]
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, Vern9(), abstol=1e-9, reltol=1e-3)
    Integral = quadgk(sol, ttr, ttr+t; rtol=1e-9, atol=0, maxevals=10^7, order=7, norm=norm)

    return Integral[1]
end

function calculate_ssquared(f, p, n, d0, miu, t=1e3, ttr=100.)
    u0 = [0, 10., 0]
    tspan = (0., ttr)
    transientprob = ODEProblem(f, u0, tspan, p)
    sol = solve(transientprob, Vern9(), abstol=1e-9, reltol=1e-3, save_everystep=false)
    tspan = (ttr, ttr+t)
    u0 = sol[end]
    prob = ODEProblem(f, u0, tspan, p)
    sol = solve(prob, Vern9(), abstol=1e-9, reltol=1e-3)
    Integral = quadgk(t->(sol(t) - miu).^2, ttr, ttr+t; rtol=1e-9, atol=0, maxevals=10^7, order=7, norm = norm)

    return Integral[1]
end

function calculate_C12(f, p, n, d0, miu, ssqr, t=1e3, ttr=100.)

    u0 = initial_conditions(f, p, n, d0)
    tspan = (ttr,ttr+t)
    prob = FunctionRepeater(f, p, u0, d0, tspan)
    function prob_func(prob, i, repeat)
        FunctionRepeater(f,p,u0[i],d0,tspan)
    end
    function output(solution, i)
        chestie = DiffEqArray([((solution[1:3,i]-miu)⋅(solution[4:6,i]-miu)) for i in 1:length(solution.t)],solution.t)
        return chestie, false
    end
    # s patrat este scalar produs scalar BOULE care esti
    monteprob = MonteCarloProblem(prob ,prob_func = prob_func, output_func = output)
    sol = solve(monteprob, Vern9(), abstol=1e-9, reltol=1e-3, parallel_type=:threads, num_monte=n, saveat= 0.01)
    intercalc = timeseries_steps_mean(sol)

    return intercalc./ssqr
end

function plot_galore1(d0, p, u0)
    prob = FunctionRepeater(f, p, u0, d0, tspan)

    function prob_func(prob, i, repeat)
    	FunctionRepeater(f,p,u0[i],d0,tspan)
    end

    monteprob = MonteCarloProblem(prob,prob_func = prob_func,
        output_func = output_func)
    sol = solve(monteprob, Vern9(), abstol=1e-9, reltol=1e-3,
        parallel_type=:threads, num_monte=n, saveat=0.001)
    dmed = timeseries_steps_mean(sol)

    return dmed
end

function plot_galore2(d0, p, u0)
    prob = FunctionRepeater(f, p, u0, d0, tspan)

    function prob_func(prob, i, repeat)
    	FunctionRepeater(f,p,u0[i],d0,tspan)
    end

    monteprob = MonteCarloProblem(prob,prob_func = prob_func,
        output_func = output_func)
    sol = solve(monteprob, Vern9(), abstol=1e-9, reltol=1e-3,
        parallel_type=:threads, num_monte=n, save_everystep = false, save_start = false)
    dmed = timeseries_steps_mean(sol)

    return dmed.u[end]
end

function plot_helper(d0, ρ)
    p = [10, ρ, 8/3]
    u0 = initial_conditions(f, p, n, d0)
    dmed = plot_galore2(d0, p, u0)
end

f = Lorenz
n = 100
p = [10, 180.7, 8/3]
# u0 = [0, 10., 0]
# u00 = [u0+rand(3) for i=1:n]
d0 = 1e-8
tspan = (0,100.)
# u0 = initial_conditions(f, p, n, d0)

#=
prob = FunctionRepeater(f, p, u0, d0, tspan)
monteprob = MonteCarloProblem(prob,prob_func = prob_func,
    output_func = output_func)
=#
# SecondLorentz = Systems.lorenz([0, 10., 0], ρ = 180.78) #Calcul lyapunov + T_λ din DynamicalSystems
# λ = lyapunov(SecondLorentz, 10000, d0 = d0)
# T_λ = (1/λ)*log(0.001/d0)

# sol = solve(monteprob, Vern9(), abstol=1e-9, reltol=1e-3, parallel_type=:threads, num_monte=n, saveat=0.001)
# dmed = timeseries_steps_mean(sol)

plt = plot()
for d0 = [1e-8, 1e-5]
    for ρ in [180.7, 181.10, 180.95]
        p = [10, ρ, 8/3]
        u0 = initial_conditions(f, p, n, d0)
        dmed = plot_galore1(d0, p, u0)
        plot!(plt, dmed.t, dmed.u, yscale=:log10)
    end
end
plt

plt2 = plot()
    for ρ in [180.7, 181.10, 180.95]
        plot!(plt2, [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4], d0->plot_helper(d0, ρ), yscale=:log10, xscale=:log10, m=4)
    end
plt2

miu = calculate_miu(f, p, n, d0)
ssqr = calculate_ssquared(f, p, n, d0, miu)
calculate_C12(f, p, n, d0, miu, ssqr)

# x = dmed.t[dmed.t .< T_λ] #comenzile pentru plotul fitarii dupa formula exponentiala
# y = dmed[dmed.t .< T_λ]
# plot(x,d0*exp.(fit_λ(sol, T_λ)*x), yscale =:log10)
#plot(x,y, legend = false, yscale=:log10, c=:blue)

#=
distantemedii = log10.(dmed.u[dmed.t .<10]) #plot fitare liniara cu linear regression (cele mai mici patrate)
a, b = linreg(timpimedii, distantemedii)
plot(timpimedii, timpimedii*b + a, c =:blue, label="")
plot!(timpimedii,distantemedii, label="$n",c=:red)
=#

# plot(dmed.t,dmed.u,yaxis=:log10)  #porcarii scrise pre implementare Monte Carlo
# sol = solve(probs,Vern9(),abstol=1e-9,reltol=1e-3)
# d = [norm(sol[1:3,i]-sol[4:6,i]) for i in 1:length(sol.t)]
# plot(sol.t,d)
