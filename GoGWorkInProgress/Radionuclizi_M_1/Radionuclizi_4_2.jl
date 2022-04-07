using Plots
using LaTeXStrings

# Graficul in timp al activitatii unei probe de Al27 
# Activata intr-un flux pulsat de neutroni
# dN/dt = ΣₐΦ - λN 

gr();

cd(@__DIR__) #Adauga path-ul la fisier

struct activitate
    Λ
    t
end

K = 2.565e7; # ΣₐΦ viteza de activare in  Hz
T_jum = 2.3 * 60; # Timp injumatatire Al28 in s
λ = log(2)/T_jum; # in Hz
τ = T_jum*2; # Intervalul fixat de timp pentru ciclurile de activare si pauza
h = 100; # Nr de intervale in care impartim timpul, elementul finit!
m = Int(100); # Nr de cicluri activare-pauza

# Calculam valorile punctuale pentru N(t) la al n-lea ciclu
function Λ_activare_n(n, t)
    Cₙ = 0.0
    for i in 0:(2*n - 2)
        Cₙ = Cₙ + (-1)^(i+1) * exp(-i*λ*τ)
    end
    return K + K * Cₙ * exp(-λ*t)
end
function Λ_pauza_n(n, t)
    Cₙ = 0.0
    for i in 0:(2*n - 1)
        Cₙ = Cₙ + (-1)^i * exp(-i*λ*τ)
    end
    return K * Cₙ * exp(-λ*t)
end

# Construim functia completa Λ(t), vector de 2 vectori
function Λ_t(m)
    Activitate = activitate(Float64[], Float64[])
    t = 0.0;
    Δt = τ/h;
    for n in 1:m
        for i in 1:h
            push!(Activitate.Λ, Λ_activare_n(n, t))
            push!(Activitate.t, t + τ * 2*(n-1))
            t = t + Δt 
        end
        t = 0.0
        for i in 1:h
            push!(Activitate.Λ, Λ_pauza_n(n, t))
            push!(Activitate.t, t + τ * (2*n-1))
            t = t + Δt 
        end    
        t = 0.0
    end
    push!(Activitate.Λ, Λ_pauza_n(m, τ))
    push!(Activitate.t, 2*m*τ)
    return Activitate
end

Λ = Λ_t(m)


plt = plot(
    Λ.t, Λ.Λ,
    xlabel = latexstring("Unități de timp în multiplii de \$\\tau\$"), 
    ylabel = L"\Lambda \: \textrm{[Bq]}",
    framestyle = :box,
    title = "Evoluția în timp a activității pentru $(m) cicluri activare-pauză",
    legend = false,
    size = (1000, 1000),
    color = :red,
    xlims = (0, maximum(Λ.t)*1.05),
    ylims = (0, maximum(Λ.Λ)*1.15)
);
for i in 1:m
    plt = plot!([2*i*τ, 2*i*τ], [0, Λ.Λ[round.(Λ.t) .== round(2*i*τ)][1]], label = "", ls = :dashdot, color = :green);
end
Δt = Int(round(2*m/10))
if Δt == 0
    Δt = 1
end
x = collect(Δt*τ:Δt*τ:τ*m*2);
y = [latexstring("\$ $(i) \\: \\tau\$") for i in Δt:Δt:2*m]; 
plt = xticks!(x, y);
plt = annotate!(Δt*τ*1.5, maximum(Λ.Λ)*1.1, latexstring("\$\\textrm{T}_{\\frac{1}{2}} = $(round(T_jum/60, digits = 2)) \\: \\textrm{[min]}\$"));
plt = annotate!(Δt*τ*1.5, maximum(Λ.Λ)*1.05, latexstring("\$\\tau \\: = $(round(τ/60, digits = 2)) \\: \\textrm{[min]}\$"));
display(plt)