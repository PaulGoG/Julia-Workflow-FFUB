using Plots
using LaTeXStrings

# Cod de calcul pentru graficele seriilor radioactive aferente temei 4

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Datele de intrare pentru cod
Λ₀ = 56.1; # μCi
T_jum_1 = 4.468*1e9; # ani
T_jum_2 = 0.06603; # ani

# Se alege functia min sau max pentru a scala graficul la timpul de injumatatire mai mic, sau mai mare dintre cei 2
T_scalare = max(T_jum_1, T_jum_2);

# Se introduc numele celor 2 radionuclizi in format LaTeX
nucleu_1 = L"^{238}\textrm{U}";
nucleu_2 = L"^{234}\textrm{Th}";

titlu = latexstring("Evoluția în timp a activității pentru seria de dezintegrări radioactive $(nucleu_1) \$\\rightarrow\$ $(nucleu_2)");

λ₁ = log(2)/T_jum_1; # ani^-1
λ₂ = log(2)/T_jum_2; # ani^-1

tₘ = (λ₂ - λ₁)^(-1) * log(λ₂/λ₁); # Valoarea t pentru care Λ isi atinge maximul

# Calculul valorilor punctuale folosind formulele deduse teoretic din ecuatiile diferentiale
function Lambda_1(t)
    return Λ₀ * exp(-λ₁*t)
end
function Lambda_2(t)
    return Λ₀ * (λ₂/(λ₂ - λ₁)) * (exp(-λ₁*t) - exp(-λ₂*t))
end

# Partea din cod care realizeaza reprezentarile grafice
x = collect(0.0:T_scalare/10:10*T_scalare); # Reglam axa timpului la 10 timpi de injumatatire
# Reprezentare Λ₁(t)
y = Lambda_1.(x);
plt = plot(
    x, 
    y,
    xlims = (minimum(x), maximum(x) * 1.01),
    ylims = (0, Λ₀ * 1.05),
    label = latexstring("$(nucleu_1):  \$\\Lambda_1(t)\$"),
    xlabel = latexstring("Unități de timp multiplii ai \$\\textrm{T}_{\\frac{1}{2}}\$, \$\\textrm{t}_{\\textrm{M}} \\simeq\$ $(round(tₘ/T_scalare, digits=2)) \$\\textrm{T}_{\\frac{1}{2}}\$"),
    ylabel = L"\Lambda \: \textrm{[\mu Ci]}",
    framestyle = :box,
    title = titlu,
    minorgrid = :true,
    size = (1000, 1000),
    legend = :topright
);
# Reprezentare Λ₂(t)
y = Lambda_2.(x);
plt = plot!(
    x, y,
    label = latexstring("$(nucleu_2):  \$\\Lambda_2(t)\$")
);
# Trasare linie verticala la valoarea tₘ
plt = plot!(
    [tₘ, tₘ], [0, Lambda_1(tₘ)],
    ls = :dashdot,
    label = "",    
);

# Scalare pe axa timpului si modificare nume diviziuni
x = collect(T_scalare:T_scalare:10*T_scalare); push!(x, tₘ);
y = [latexstring("\$\\textrm{T}_{\\frac{1}{2}}\$")];
append!(y, [latexstring("\$ $(i) \\: \\textrm{T}_{\\frac{1}{2}}\$") for i in 2:10]); 
push!(y, L"t_M");
plt = xticks!(x, y);

# Legenda reprezentarii grafice
plt = annotate!(5*T_scalare, Λ₀*1, latexstring("\$\\textrm{T}_{\\frac{1}{2}}\$ $(nucleu_1) = $(T_jum_1) ani"));
plt = annotate!(5*T_scalare, Λ₀*0.95, latexstring("\$\\textrm{T}_{\\frac{1}{2}}\$ $(nucleu_2) = $(T_jum_2) ani"));
plt = annotate!(5*T_scalare, Λ₀*0.9, latexstring("\$\\textrm{t}_{M}\$ = \$$(round(tₘ, digits=2))\$ ani"));

display(plt);
#savefig(plt, "Grafice\\Evolutie_Activitate_$(nucleu_1)_$(nucleu_2).png")
#savefig(plt, "Grafice/Evolutie_Activitate_$(nucleu_1)_$(nucleu_2).png")