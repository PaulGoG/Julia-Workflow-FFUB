using Distributed
addprocs(4)
@everywhere using DifferentialEquations
using Makie

#Differential Equations of Motion for one Ball
@everywhere function f(du,u,p,t)
  du[1] = u[2]
  du[2] = -p
  du[3] = u[4]
  du[4] = 0.0
  du[5] = u[6]
  du[6] = 0.0
end

#Conditon of collision with the walls for one Ball
@everywhere function condition(out,u,t,integrator)
  out[1] = (u[1] - 9.5)*(u[1]-0.5)
  out[2] = (u[3] - 9.5)*(u[3]-0.5)
  out[3] = (u[5] - 9.5)*(u[5]-0.5)
end

#Reflection for one ball
@everywhere function affect!(integrator, idx)
  if idx == 1
    integrator.u[2] = -0.9*integrator.u[2]
  elseif idx == 2
    integrator.u[4] = -0.9integrator.u[4]
  elseif idx == 3
    integrator.u[6] = -0.9integrator.u[6]
  end
end

#Callback vector
cb = VectorContinuousCallback(condition,affect!,3)
#No. of balls to simulate
N = 10

u0 = [5.0,10*(2*rand()-1),5.0,10*(2*rand()-1),5.0,10*(2*rand()-1)]
tspan = (0.0,10.0)
p = 9.8

prob = ODEProblem(f,u0,tspan,p)

@everywhere function prob_func(prob,i,repeat)
  ODEProblem(prob.f,[5.0,10*(2*rand()-1),5.0,10*(2*rand()-1),5.0,10*(2*rand()-1)],prob.tspan,prob.p)
end

ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
@time sim = solve(ensemble_prob,Vern9(),EnsembleDistributed(),callback=cb,dt=1e-2,adaptive=false,trajectories=N)

#Animation part
t = Node(0.)

ts = tspan[1]:0.01:tspan[end]

p = lift(t; init = [Point3f0(sim[1](t[])[1],sim[1](t[])[3],sim[1](t[])[5])]) do val
  p[] = [Point3f0(sim[i](val)[3],sim[i](val)[5],sim[i](val)[1]) for i in 1:N]
end

limits = FRect3D((0,0,0),(10,10,10))

scene = scatter(p[],limits=limits,markersize=0.5)

for i in ts
  push!(t,i)
  sleep(1/40)
end

#to record the video
record(scene,"animated_ball.mp4",range(tspan[1],stop=tspan[end],step=0.01)) do i
  push!(t,i)
  sleep(1/120)
end
