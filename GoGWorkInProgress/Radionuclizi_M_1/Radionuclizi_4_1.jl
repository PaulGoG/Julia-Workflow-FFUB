using Plots
using LaTeXStrings

# Graficul activitatilor in timp al primilor doi izotopi radioactivi 
# Dintr-un lant de dezintegrari succesive 1 -> 2 -> 3 ...

gr();

cd(@__DIR__) #Adauga path-ul la fisier

Λ₀ = 56.1; # μCi
T_jum_1 = 4.468 * 10^9; # ani
T_jum_2 = 0.066; # ani
T_scalare = max(T_jum_1, T_jum_2); # se alege functia min sau max in functie de scalarea graficului

nucleu_1 = L"^{238}\textrm{U}";
nucleu_2 = L"^{234}\textrm{Th}";
titlu = latexstring("Graficul evoluției în timp a activității pentru seria de dezintegrări radioactive $(nucleu_1) \$\\rightarrow\$ $(nucleu_2)");

λ₁ = log(2)/T_jum_1; # ani^-1
λ₂ = log(2)/T_jum_2; # ani^-1

tₘ = (λ₂ - λ₁)^(-1) * log(λ₂/λ₁);

function Lambda_1(t)
    return Λ₀ * exp(-λ₁*t)
end
function Lambda_2(t)
    return Λ₀ * (λ₂/(λ₂ - λ₁)) * (exp(-λ₁*t) - exp(-λ₂*t))
end

x = collect(0.0:T_scalare/10:10*T_scalare);
y = Lambda_1.(x);
plt = plot(
    x, y,
    xlims = (minimum(x), maximum(x) * 1.01),
    ylims = (0, Λ₀ * 1.05),
    label = L"\Lambda_1(t)",
    xlabel = latexstring("Unități de timp multiplii ai \$\\textrm{T}_{\\frac{1}{2}}\$, t\$_M\$ = $(round(tₘ, digits=2)) ani, sau $(round(tₘ/T_scalare, digits=2)) \$\\textrm{T}_{\\frac{1}{2}}\$"),
    ylabel = L"\Lambda \: \textrm{[\mu Ci]}",
    framestyle = :box,
    title = titlu,
    minorgrid = :true,
    size = (1000, 1000),
    legend = :bottomleft
);
y = y = Lambda_2.(x);
plt = plot!(
    x, y,
    label = L"\Lambda_2(t)"
);
plt = plot!(
    [tₘ, tₘ], [0, Lambda_1(tₘ)],
    ls = :dashdot,
    label = "",    
    );

x = collect(T_scalare:T_scalare:10*T_scalare); push!(x, tₘ);
y = [latexstring("\$\\textrm{T}_{\\frac{1}{2}}\$")];
append!(y, [latexstring("\$ $(i) \\: T_{\\frac{1}{2}}\$") for i in 2:10]); 
push!(y, L"t_M");
plt = xticks!(x, y);

plt = annotate!(5*T_scalare, Lambda_1(tₘ)*0.9, latexstring("\$\\textrm{T}_{\\frac{1}{2}}\$ $(nucleu_1) = $(T_jum_1) ani"));
plt = annotate!(5*T_scalare, Lambda_1(tₘ)*0.85, latexstring("\$\\textrm{T}_{\\frac{1}{2}}\$ $(nucleu_2) = $(T_jum_2) ani"));

display(plt);
# savefig(plt, "Grafice\\Evolutie_Activitate_$(nucleu_1)_$(nucleu_2).png")