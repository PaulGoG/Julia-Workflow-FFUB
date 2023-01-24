using Plots
using LaTeXStrings

# Cod de calcul pentru compunerea a doua rezonante Breit-Wigner avand aceleasi numere cuantice (coerenta)

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

struct distributie
    E
    B_1
    B_2
    P_1
    P_2
end

function Bₖ(E, Eₖ, Γₖ)
    return (0.5 * Γₖ)/(Eₖ - E - 0.5*im*Γₖ)
end
function dσ_dε(B₁, B₂, ϕ)
    return abs2(B₁ + B₂*exp(im*ϕ))*0.5
end

E = collect(2000.0:0.1:3000.0);
P_mass = distributie(Float64[], Float64[], Float64[], Float64[], Float64[]);
ε₁ = 2300;
Γ₁ = 150;
ε₂ = 2340;
Γ₂ = 320;
ϕ₃₀ = π/6;
ϕ₄₅ = π/4;


for i in eachindex(E)
    B₁ = Bₖ(E[i], ε₁, Γ₁)
    B₂ = Bₖ(E[i], ε₂, Γ₂)
    push!(P_mass.E, E[i])
    push!(P_mass.B_1, abs2(B₁))
    push!(P_mass.B_2, abs2(B₂))
end
for i in eachindex(E)
    B₁ = Bₖ(E[i], ε₁, Γ₁)
    B₂ = Bₖ(E[i], ε₂, Γ₂)    
    push!(P_mass.P_1, dσ_dε(B₁, B₂, ϕ₃₀))
    push!(P_mass.P_2, dσ_dε(B₁, B₂, ϕ₄₅))
end
plot(P_mass.E, P_mass.B_1);
plot!(P_mass.E, P_mass.B_2);
plot!(P_mass.E, P_mass.P_1);
plot!(P_mass.E, P_mass.P_2)