using Plots
using Trapz
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
    return abs2(B₁ + B₂*exp(im*ϕ))
end
function Grafic(x, y, titlu, eticheta, culoare, axa_x, axa_y)
    plt = plot(
        x,
        y,
        size = (1600, 1200),
        minorgrid = true,
        framestyle = :box,
        xlims = (minimum(x), maximum(x)),
        yformatter = :scientific,
        xlabel = axa_x,
        ylabel = axa_y,
        title = titlu,
        label = eticheta,
        color = culoare
    )
    indice_min = findfirst(i -> minimum(y) >= abs(maximum(y)/2 - i), y)
    indice_max = findlast(i -> minimum(y) >= abs(maximum(y)/2 - i), y)
    x = [x[indice_min], x[indice_max]]
    y = [y[indice_min], y[indice_max]]
    plot!(
        plt,
        x,
        y,
        ls = :dashdot,
        color = culoare,
        label = ""
    )
    return plt
end
function Grafic(x, y, eticheta, culoare, plt)
    plot!(
        plt,
        x,
        y,
        label = eticheta,
        color = culoare
    )
    indice_min = findfirst(i -> minimum(y) >= abs(maximum(y)/2 - i), y)
    indice_max = findlast(i -> minimum(y) >= abs(maximum(y)/2 - i), y)
    x = [x[indice_min], x[indice_max]]
    y = [y[indice_min], y[indice_max]]
    plot!(
        plt,
        x,
        y,
        ls = :dashdot,
        color = culoare,
        label = ""
    )
    return plt
end
function Afisare_grafic(plt, titlu)
    display(plt)
    #savefig("$titlu.png")
end

E = 1600.0:0.1:3000.0;
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
Plot_B_nenormat = Grafic(P_mass.E, P_mass.B_1, "Distributiile Breit-Wigner nenormate", L"\mathrm{B}_1", :blue, "E [MeV]", "P(E)");
Plot_B_nenormat = Grafic(P_mass.E, P_mass.B_2,L"\mathrm{B}_2", :red, Plot_B_nenormat);
C_1 = sqrt(1/trapz(P_mass.E, P_mass.B_1));
C_2 = sqrt(1/trapz(P_mass.E, P_mass.B_2));
P_mass.B_1 .*= C_1^2;
P_mass.B_2 .*= C_2^2;
Plot_B_normat = Grafic(P_mass.E, P_mass.B_1, "Distributiile Breit-Wigner normate", L"\mathrm{B}_1", :blue, "E [MeV]", "P(E)");
Plot_B_normat = Grafic(P_mass.E, P_mass.B_2, L"\mathrm{B}_2", :red, Plot_B_normat);
for i in eachindex(E)
    B₁ = Bₖ(E[i], ε₁, Γ₁) * C_1
    B₂ = Bₖ(E[i], ε₂, Γ₂) * C_2
    push!(P_mass.P_1, dσ_dε(B₁, B₂, ϕ₃₀))
    push!(P_mass.P_2, dσ_dε(B₁, B₂, ϕ₄₅))
end
C_30 = sqrt(1/trapz(P_mass.E, P_mass.P_1));
C_45 = sqrt(1/trapz(P_mass.E, P_mass.P_2));
P_mass.P_1 .*= C_30^2;
P_mass.P_2 .*= C_45^2;

Plot_B_P_30 = Grafic(
    P_mass.E, 
    P_mass.P_1, 
    latexstring("Distributiile Breit-Wigner normate & distributia compusa coerent cu un unghi de defazaj  \$\\varphi = 30^o\$"), 
    L"\left|\mathrm{C}_1 \mathrm{B}_1 + \mathrm{C}_2 \mathrm{B}_2 \mathrm{e}^{i \varphi}\right|^2",
    :black, "E [MeV]", "P(E)"
);
Plot_B_P_30 = Grafic(P_mass.E, P_mass.B_1, L"\mathrm{B}_1", :blue, Plot_B_P_30);
Plot_B_P_30 = Grafic(P_mass.E, P_mass.B_2, L"\mathrm{B}_2", :red, Plot_B_P_30);

Plot_B_P_45 = Grafic(
    P_mass.E, 
    P_mass.P_2, 
    latexstring("Distributiile Breit-Wigner normate & distributia compusa coerent cu un unghi de defazaj  \$\\varphi = 45^o\$"), 
    L"\left|\mathrm{C}_1 \mathrm{B}_1 + \mathrm{C}_2 \mathrm{B}_2 \mathrm{e}^{i \varphi}\right|^2",
    :black, "E [MeV]", "P(E)"
);
Plot_B_P_45 = Grafic(P_mass.E, P_mass.B_1, L"\mathrm{B}_1", :blue, Plot_B_P_45);
Plot_B_P_45 = Grafic(P_mass.E, P_mass.B_2, L"\mathrm{B}_2", :red, Plot_B_P_45);

Plot_P_30_45 = Grafic(
    P_mass.E, 
    P_mass.P_1, 
    latexstring("Distributiile Breit-Wigner compuse coerent cu unghiurile de defazaj  \$\\varphi_1 = 30^o\$ & \$\\varphi_2 = 45^o\$"), 
    L"\varphi_1 = 30^o",
    :green, "E [MeV]", "P(E)"
);
Plot_P_30_45 = Grafic(P_mass.E, P_mass.P_2, L"\varphi_2 = 45^o", :orange, Plot_P_30_45);

Afisare_grafic(Plot_B_nenormat, "B_nenormat");
Afisare_grafic(Plot_B_normat, "B_normat");
Afisare_grafic(Plot_B_P_30, "P_30_B");
Afisare_grafic(Plot_B_P_45, "P_45_B");
Afisare_grafic(Plot_P_30_45, "P_30_45");