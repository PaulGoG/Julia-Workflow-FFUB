# Monte Carlo simulation of a simple pendulum's phase space

using Plots
using DifferentialEquations

function Pendulum(du, u, p, t)
    du[1] = u[2]
    du[2] = -sin(u[1])
end

function RepeatFunction(f, u0, tspan)
    ODEProblem(f, u0, tspan)
end

function prob_func(prob, i, repeat)
    u0 = [rand()*2*π - π, rand()]
    RepeatFunction(f, u0, tspan)
end

f = Pendulum
tspan = (0,100.)
u0 = [rand()*2*π - π, rand()]
n = 100

prob = ODEProblem(f, u0, tspan)
monteprob = MonteCarloProblem(prob, prob_func = prob_func)
sol = solve(monteprob, SSPRK432(), parallel_type=:threads,
 num_monte = n)

plot(sol, vars = (1,2), xlims = (-π, π) )

# sol = solve(prob, SSPRK432())
# plot(sol, vars = (2,1))
