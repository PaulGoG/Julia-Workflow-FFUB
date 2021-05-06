#=
    i∂ₜΨ = -1/2 ∂ₓₓΨ - |Ψ|²Ψ by MOL
    n -> indice temporal
    j -> indice spatial
    conditii de frontiera spatiale periodice
=#

using Plots; #plotlyjs();
using Distributions

function Sⁿ(Ψⁿⱼ₋₁, Ψⁿⱼ, Ψⁿⱼ₊₁, dx)
    return (im/2) * ((Ψⁿⱼ₊₁ - 2*Ψⁿⱼ + Ψⁿⱼ₋₁)/(dx^2)) + im*Ψⁿⱼ*abs2(Ψⁿⱼ)
end

function RungeKutta3(Ψⁿⱼ₋₁, Ψⁿⱼ, Ψⁿⱼ₊₁, dx, dt)

    k1 = Sⁿ(Ψⁿⱼ₋₁, Ψⁿⱼ, Ψⁿⱼ₊₁, dx)
    k2 = Sⁿ(Ψⁿⱼ₋₁ + dt*k1/2, Ψⁿⱼ + dt*k1/2, Ψⁿⱼ₊₁ + dt*k1/2, dx)
    k3 = Sⁿ(Ψⁿⱼ₋₁ + dt*k2/2, Ψⁿⱼ + dt*k2/2, Ψⁿⱼ₊₁ + dt*k2/2, dx)

    return Ψⁿⱼ + dt/6 * (k1 + k2*2 + k3)
end

a = 10;
timp = 200;
#=
dx = 2*a/(N_x-1); # x de la -a la a
dt = t/(N_timp-1); # t de la 0 la t
N_timp = Int(1e4);
N_x = Int(1e2);
=#

dx = 0.5;
dt = 10^(-2) * dx^2;
x = collect(-a:dx:a);
t = collect(0:dt:timp);
N_timp = length(t);
N_x = length(x);
Ψ = zeros(N_timp, N_x);
Ψ = ComplexF64.(Ψ);
#Ψ[1, :] = [pdf(Normal(0, 1), x[i]) for i in 1:length(x)]; conditie initiala Gaussiana
#Cond. initiala ecuatia analitica din articol
for i in 1:Int(floor(length(x)/2))
    Ψ[1, i] = 2 * exp(-im*0.1*x[i])/cosh(2*(x[i] + 5))
end
for i in Int(floor(length(x)/2)):length(x)
    Ψ[1, i] = 2 * exp(im*0.1*x[i])/cosh(2*(x[i] - 5))
end


for n in 1:(N_timp-1)
    Ψ[n+1, 1] = RungeKutta3(Ψ[n, N_x - 1], Ψ[n, 1], Ψ[n, 2], dx, dt)
    Ψ[n+1, N_x] = RungeKutta3(Ψ[n, N_x - 1], Ψ[n, N_x], Ψ[n, 1], dx, dt)
    for j in 2:(N_x - 1)
        Ψ[n+1, j] = RungeKutta3(Ψ[n, j-1], Ψ[n, j], Ψ[n,j+1], dx, dt)
    end
end

ρ = abs2.(Ψ)
surface(x, t, ρ,
xlabel = "x", ylabel = "t",
ylims = (0, t[N_timp]),
xlims = (-7.5, 7.5))