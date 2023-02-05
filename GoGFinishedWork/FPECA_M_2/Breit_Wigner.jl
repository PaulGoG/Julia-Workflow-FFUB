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
        size = (1000, 1000),
        dpi = 600,
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
    indice_min = findfirst(i -> minimum(y)*0.005 >= abs(maximum(y)/2 - i), y)
    indice_max = findlast(i -> minimum(y)*0.005 >= abs(maximum(y)/2 - i), y)
    x = [x[indice_min], x[indice_max]]
    y = [y[indice_min], y[indice_max]]
    plot!(
        plt,
        x,
        y,
        ls = :dot,
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
    indice_min = findfirst(i -> minimum(y)*0.005 >= abs(maximum(y)/2 - i), y)
    indice_max = findlast(i -> minimum(y)*0.005 >= abs(maximum(y)/2 - i), y)
    x = [x[indice_min], x[indice_max]]
    y = [y[indice_min], y[indice_max]]
    plot!(
        plt,
        x,
        y,
        ls = :dot,
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
Plot_B_nenormat = Grafic(P_mass.E, P_mass.B_1, "Distribuțiile Breit-Wigner nenormate", L"|\mathrm{B}_1|^2(\mathrm{E}_1 = 2300 \: & \: \Gamma_1 = 150)", :blue, "E [MeV]", "P(E)");
Plot_B_nenormat = Grafic(P_mass.E, P_mass.B_2, L"|\mathrm{B}_2|^2(\mathrm{E}_2 = 2340 \: & \: \Gamma_2 = 320)", :orange, Plot_B_nenormat);
Afisare_grafic(Plot_B_nenormat, "B_nenormat");
C_1 = sqrt(1/trapz(P_mass.E, P_mass.B_1));
C_2 = sqrt(1/trapz(P_mass.E, P_mass.B_2));
P_mass.B_1 .*= C_1^2;
P_mass.B_2 .*= C_2^2;
Plot_B_normat = Grafic(P_mass.E, P_mass.B_1, "Distribuțiile Breit-Wigner normate", L"|\mathrm{C}_1\mathrm{B}_1|^2(\mathrm{E}_1 = 2300 \: & \: \Gamma_1 = 150)", :blue, "E [MeV]", "P(E)");
Plot_B_normat = Grafic(P_mass.E, P_mass.B_2, L"|\mathrm{C}_2\mathrm{B}_2|^2(\mathrm{E}_2 = 2340 \: & \: \Gamma_2 = 320)", :orange, Plot_B_normat);
Afisare_grafic(Plot_B_normat, "B_normat");
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
    latexstring("Distribuțiile Breit-Wigner normate ce interferă coerent cu defazajul  \$\\varphi = 30^o\$"), 
    L"|\mathrm{C}|^2 \; |\mathrm{C}_1 \mathrm{B}_1 + \mathrm{C}_2 \mathrm{B}_2 \mathrm{e}^{i \varphi}|^2",
    :black, "E [MeV]", "P(E)"
);
Plot_B_P_30 = Grafic(P_mass.E, P_mass.B_1, L"|\mathrm{C}_1\mathrm{B}_1|^2(\mathrm{E}_1 = 2300 \: & \: \Gamma_1 = 150)", :orange, Plot_B_P_30);
Plot_B_P_30 = Grafic(P_mass.E, P_mass.B_2, L"|\mathrm{C}_2\mathrm{B}_2|^2(\mathrm{E}_2 = 2340 \: & \: \Gamma_2 = 320)", :yellow, Plot_B_P_30);
Afisare_grafic(Plot_B_P_30, "P_30_B");

Plot_B_P_45 = Grafic(
    P_mass.E, 
    P_mass.P_2, 
    latexstring("Distribuțiile Breit-Wigner normate ce interferă coerent cu defazajul  \$\\varphi = 45^o\$"), 
    L"|\mathrm{C}|^2 \; |\mathrm{C}_1 \mathrm{B}_1 + \mathrm{C}_2 \mathrm{B}_2 \mathrm{e}^{i \varphi}|^2",
    :black, "E [MeV]", "P(E)"
);
Plot_B_P_45 = Grafic(P_mass.E, P_mass.B_1, L"|\mathrm{C}_1\mathrm{B}_1|^2(\mathrm{E}_1 = 2300 \: & \: \Gamma_1 = 150)", :orange, Plot_B_P_45);
Plot_B_P_45 = Grafic(P_mass.E, P_mass.B_2, L"|\mathrm{C}_2\mathrm{B}_2|^2(\mathrm{E}_2 = 2340 \: & \: \Gamma_2 = 320)", :yellow, Plot_B_P_45);
Afisare_grafic(Plot_B_P_45, "P_45_B");

indice_min = findfirst(i -> minimum(P_mass.P_1)*0.005 >= abs(maximum(P_mass.P_1)/2 - i), P_mass.P_1);
indice_max = findlast(i -> minimum(P_mass.P_1)*0.005 >= abs(maximum(P_mass.P_1)/2 - i), P_mass.P_1);
Γ_P_1 = round(P_mass.E[P_mass.P_1 .== P_mass.P_1[indice_max]][1] - P_mass.E[P_mass.P_1 .== P_mass.P_1[indice_min]][1], digits = 2);
E_P_1 = round(P_mass.E[P_mass.P_1 .== maximum(P_mass.P_1)][1], digits = 2);
Plot_P_30_45 = Grafic(
    P_mass.E, 
    P_mass.P_1, 
    latexstring("Distribuțiile Breit-Wigner de interferență la unghiurile de defazaj  \$\\varphi_a = 30^o\$ & \$\\varphi_b = 45^o\$"), 
    latexstring("\$\\varphi_a = 30^o\$, \$\\mathrm{E}_a = $(E_P_1)\$, \$\\Gamma_a = $(Γ_P_1)\$"),
    :black, "E [MeV]", "P(E)"
);
indice_min = findfirst(i -> minimum(P_mass.P_2)*0.005 >= abs(maximum(P_mass.P_2)/2 - i), P_mass.P_2);
indice_max = findlast(i -> minimum(P_mass.P_2)*0.005 >= abs(maximum(P_mass.P_2)/2 - i), P_mass.P_2);
Γ_P_2 = round(P_mass.E[P_mass.P_2 .== P_mass.P_2[indice_max]][1] - P_mass.E[P_mass.P_2 .== P_mass.P_2[indice_min]][1], digits = 2);
E_P_2 = round(P_mass.E[P_mass.P_2 .== maximum(P_mass.P_2)][1], digits = 2);
Plot_P_30_45 = Grafic(P_mass.E, P_mass.P_2, latexstring("\$\\varphi_b = 45^o\$, \$\\mathrm{E}_b = $(E_P_2)\$, \$\\Gamma_b = $(Γ_P_2)\$"), :orange, Plot_P_30_45);
Afisare_grafic(Plot_P_30_45, "P_30_45");