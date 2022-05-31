using Plots
using LaTeXStrings

# Cod de calcul pentru graficul in timp al activitatii unei probe de Al-27 
# Activata intr-un flux pulsat de neutroni conform ecuatiei
# dN/dt = ΣₐΦ - λN, aferent temei 4 

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

struct Activitate
    Λ
    t
end

# Date de intrare pentru cod
K = 2.565e7; # ΣₐΦ viteza de activare in  s^-1
T_jum = 2.3 * 60; # Timp injumatatire Al-28 in s

# Intervalul de timp pentru ciclurile de activare si pauza
τ = T_jum*5; 
# Nr de cicluri activare-pauza simulate
m = 3;
m = Int(m);
# Nr de intervale in care impartim timpul, elementul finit!
h = 100;

λ = log(2)/T_jum; # in s^-1

# Calculam valorile punctuale pentru N(t) la al n-lea ciclu folosind formulele deduse analitic
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

# Construim functia completa Λ(t) din valorile calculate punctual
function Λ_t(m)
    activitate = Activitate(Float64[], Float64[])
    t = 0.0;
    Δt = τ/h;
    for n in 1:m
        for i in 1:h
            push!(activitate.Λ, Λ_activare_n(n, t))
            push!(activitate.t, t + τ * 2*(n-1))
            t = t + Δt 
        end
        t = 0.0;
        for i in 1:h
            push!(activitate.Λ, Λ_pauza_n(n, t))
            push!(activitate.t, t + τ * (2*n-1))
            t = t + Δt 
        end    
        t = 0.0;
    end
    push!(activitate.Λ, Λ_pauza_n(m, τ))
    push!(activitate.t, 2*m*τ)
    return activitate
end

# Apelarea functiei de constructie a Λ(t)
Λ = Λ_t(m);

# Constructia reprezentarilor grafice
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
# Delimitarea pe grafic cu linie punctata a ciclurilor activare-pauza
for i in 1:m
    plt = plot!([2*i*τ, 2*i*τ], [0, Λ.Λ[round.(Λ.t) .== round(2*i*τ)][1]], label = "", ls = :dashdot, color = :green);
end
# Tratare caz numeric limita m = 1
Δt = Int(round(2*m/10))
if Δt == 0
    Δt = 1
end

# Scalare pe axa timpului si modificare nume diviziuni
x = collect(Δt*τ:Δt*τ:τ*m*2);
y = [latexstring("\$ $(i) \\: \\tau\$") for i in Δt:Δt:2*m]; 
plt = xticks!(x, y);

# Legenda reprezentarii grafice
plt = annotate!(Δt*τ*1.5, maximum(Λ.Λ)*1.1, latexstring("\$\\textrm{T}_{\\frac{1}{2}} = $(round(T_jum/60, digits = 2)) \\: \\textrm{[min]}\$"));
plt = annotate!(Δt*τ*1.5, maximum(Λ.Λ)*1.05, latexstring("\$\\tau \\: = $(round(τ/60, digits = 2)) \\: \\textrm{[min]}\$"));

display(plt)
#savefig(plt, "Grafice\\Cicluri_Activare_Pauza.png")
#savefig(plt, "Grafice/Cicluri_Activare_Pauza.png")