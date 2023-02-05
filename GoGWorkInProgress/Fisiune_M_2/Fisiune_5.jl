using Plots
using CSV
using DataFrames
using LaTeXStrings
using QuadGK
using Optim
using Trapz

# Cod de calcul pentru fitul datelor experimentale ale unor spectre de neutroni cu o distributie Maxwell

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Citire fisiere de date
dn_Gook_SL = CSV.read("Data_files/Date_experimentale/Spectru_n/U5SPGOOK.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["E", "N", "σN"], skipto = 2);
dn_Vorobyev_SL = CSV.read("Data_files/Date_experimentale/Spectru_n/U5SPVORO.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["E", "N", "σN"], skipto = 2);
dn_Gook_SCM_LF = CSV.read("Data_files/Date_experimentale/Spectru_n/U5SPCMLF.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["E", "N", "σN"], skipto = 2);
dn_Gook_SCM_HF = CSV.read("Data_files/Date_experimentale/Spectru_n/U5SPCMHF.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["E", "N", "σN"], skipto = 2);
#####
# Nₘ(E, Tₘ)
function N_Maxwell(E, T_M)
    return (2/sqrt(π)) * T_M^(-3/2) * sqrt(E) * exp(-E/T_M)
end
# Aria de sub graficul datelor experimentale pentru domeniul acoperit de punctele experimentale
function Arie_date_experimentale(dn)
    return trapz(dn.E, dn.N)
end
# Aria de sub grafic pentru domeniul acoperit de punctele experimentale a distributiei Maxwell
function Arie_teoretica(E_min, E_max, T_M)
    Maxwell(E) = (2/sqrt(π)) * T_M^(-3/2) * sqrt(E) * exp(-E/T_M)
    return quadgk(Maxwell, E_min, E_max)[1]
end
# Factorul de normare f = A_th/A_exp
function Normare_Maxwell(dn, T_M)
    E_min = first(dn.E)
    E_max = last(dn.E)
    return Arie_teoretica(E_min, E_max, T_M)/Arie_date_experimentale(dn)
end
# Optimizarea parametrica a lui Tₘ si calculul χ²
function Minimizare_Hi_patrat(dn)
    χ²(T_M) = (1/length(dn.E)) * sum([((dn.N[i] - (N_Maxwell(dn.E[i], T_M)/Normare_Maxwell(dn, T_M)))/dn.σN[i])^2 for i in eachindex(dn.N)]);
    result = optimize(χ², 0.0, 5.0);
    return [Optim.minimizer(result), Optim.minimum(result)]
end
# Normarea datelor experimentale la spectrul teoretic
function Normare_Spectru(dn, T_M)
    f = Normare_Maxwell(dn, T_M)
    for i in eachindex(dn.E)
        dn.N[i] *= f/N_Maxwell(dn.E[i], T_M)
        dn.σN[i] *= f/N_Maxwell(dn.E[i], T_M)
    end
    return dn
end
# Aici se opreste partea de calcul a programului
#####
# Constructia reprezentarilor grafice
function Grafic_plot(x, y, titlu, eticheta, axa_x, axa_y, scalare_inf, scalare_sup, scala_x, scala_y, x_maxim, culoare)
    plt = plot(
        x, 
        y, 
        xlims = (minimum(x), x_maxim),
        ylims = (minimum(y)*scalare_inf, maximum(y)*scalare_sup),
        xlabel = axa_x, 
        ylabel = axa_y, 
        framestyle = :box,
        label = eticheta,
        title = titlu,
        xscale = scala_x,
        yscale = scala_y,
        minorgrid = :true,
        size = (1000, 1000),
        dpi = 600,
        color = culoare
    )
    return plt
end
function Grafic_scatter(x, y, σ, culoare, eticheta, plt)
    scatter!(plt,
        x,
        y,
        yerr = σ,
        color = culoare,
        label = eticheta,
    )
end
function Grafic_textbox(x, y, text, plt)
    annotate!(plt,
    x,
    y,
    text
    )
    return plt
end
function Grafic_afisare(plt, titlu)
    display(plt)
    savefig(plt, "Grafice/T5_$(titlu).png")
end
#####
# Apelarea functiilor definite pentru executia programului
E = 1e-1:1e-2:20;

dn = copy(dn_Gook_SL);
T_M, χ_sq = Minimizare_Hi_patrat(dn);
Plot = Grafic_plot(E, N_Maxwell.(E, T_M), "Fitul datelor experimentale ale spectrului neutronilor prompți cu o distribuție Maxwell", latexstring("Distribuția Maxwell cu \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), "E [MeV]", "N(E)", 1.0, 2.0, :identity, :log10, maximum(dn.E)*1.05, :red);
Plot = Grafic_scatter(dn.E, dn.N, dn.σN, :black, "Date experimentale ale lui Gook în SL", Plot);
Plot = Grafic_textbox(10.25, 0.18, latexstring("\$\\chi^2\$ = $(round(χ_sq, digits = 3)), \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), Plot);
yticks!(Plot, [10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1]);
Grafic_afisare(Plot, "T_5_Gook_SL_raw");
dn_normat = Normare_Spectru(copy(dn), T_M);
Plot = Grafic_plot(E, [1.0 for i in eachindex(E)], "Fitul datelor experimentale ale spectrului neutronilor prompți cu o distribuție Maxwell", "", "E [MeV]", "Raportul spectrului neutronilor prompți la distribuția Maxwell", 0.0, 1.1, :log10, :identity, maximum(dn.E)*1.05, :red);
Plot = Grafic_scatter(dn_normat.E, dn_normat.N, dn_normat.σN, :black, "Date experimentale ale lui Gook în SL", Plot);
Plot = Grafic_textbox(4.7, 0.09, latexstring("\$\\chi^2\$ = $(round(χ_sq, digits = 3)), \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), Plot);
Grafic_afisare(Plot, "T_5_Gook_SL_normat");

dn = copy(dn_Vorobyev_SL);
T_M, χ_sq = Minimizare_Hi_patrat(dn);
Plot = Grafic_plot(E, N_Maxwell.(E, T_M), "Fitul datelor experimentale ale spectrului neutronilor prompți cu o distribuție Maxwell", latexstring("Distribuția Maxwell cu \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), "E [MeV]", "N(E)", 1.0, 3.0, :identity, :log10, maximum(dn.E)*1.05, :red);
Plot = Grafic_scatter(dn.E, dn.N, dn.σN, :black, "Date experimentale ale lui Vorobyev în SL", Plot);
Plot = Grafic_textbox(13.3, 0.28, latexstring("\$\\chi^2\$ = $(round(χ_sq, digits = 3)), \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), Plot);
yticks!(Plot, [10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1]);
Grafic_afisare(Plot, "T_5_Voro_SL_raw");
dn_normat = Normare_Spectru(copy(dn), T_M);
Plot = Grafic_plot(E, [1.0 for i in eachindex(E)], "Fitul datelor experimentale ale spectrului neutronilor prompți cu o distribuție Maxwell", "", "E [MeV]", "Raportul spectrului neutronilor prompți la distribuția Maxwell", 0.0, 1.2, :log10, :identity, maximum(dn.E)*1.05, :red);
Plot = Grafic_scatter(dn_normat.E, dn_normat.N, dn_normat.σN, :black, "Date experimentale ale lui Vorobyev în SL", Plot);
Plot = Grafic_textbox(2.8e-1, 0.1, latexstring("\$\\chi^2\$ = $(round(χ_sq, digits = 3)), \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), Plot);
Grafic_afisare(Plot, "T_5_Voro_SL_normat");

dn = copy(dn_Gook_SCM_LF);
T_M, χ_sq = Minimizare_Hi_patrat(dn);
Plot = Grafic_plot(E, N_Maxwell.(E, T_M), "Fitul datelor experimentale ale spectrului neutronilor prompți cu o distribuție Maxwell", latexstring("Distribuția Maxwell cu \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), "E [MeV]", "N(E)", 1e6, 2, :identity, :log10, maximum(dn.E)*1.05, :red);
Plot = Grafic_scatter(dn.E, dn.N, dn.σN, :black, "Date experimentale ale lui Gook în SCM pentru LF", Plot);
Plot = Grafic_textbox(5.05, 0.5, latexstring("\$\\chi^2\$ = $(round(χ_sq, digits = 3)), \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), Plot);
yticks!(Plot, [10^-4, 10^-3, 10^-2, 10^-1, 1]);
Grafic_afisare(Plot, "T_5_Gook_SCM_LF_raw");
dn_normat = Normare_Spectru(copy(dn), T_M);
Plot = Grafic_plot(E, [1.0 for i in eachindex(E)], "Fitul datelor experimentale ale spectrului neutronilor prompți cu o distribuție Maxwell", "", "E [MeV]", "Raportul spectrului neutronilor prompți la distribuția Maxwell", 0.0, 4.0, :log10, :identity, 10.0, :red);
Plot = Grafic_scatter(dn_normat.E, dn_normat.N, dn_normat.σN, :black, "Date experimentale ale lui Gook în SCM pentru LF", Plot);
Plot = Grafic_textbox(2.4, 3.65, latexstring("\$\\chi^2\$ = $(round(χ_sq, digits = 3)), \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), Plot);
Grafic_afisare(Plot, "T_5_Gook_SCM_LF_normat");

dn = copy(dn_Gook_SCM_HF);
T_M, χ_sq = Minimizare_Hi_patrat(dn);
Plot = Grafic_plot(E, N_Maxwell.(E, T_M), "Fitul datelor experimentale ale spectrului neutronilor prompți cu o distribuție Maxwell", latexstring("Distribuția Maxwell cu \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), "E [MeV]", "N(E)", 5e5, 1.5, :identity, :log10, maximum(dn.E)*1.05, :red);
Plot = Grafic_scatter(dn.E, dn.N, dn.σN, :black, "Date experimentale ale lui Gook în SCM pentru HF", Plot);
Plot = Grafic_textbox(5.18, 0.38, latexstring("\$\\chi^2\$ = $(round(χ_sq, digits = 3)), \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), Plot);
yticks!(Plot, [10^-3, 10^-2, 10^-1, 1]);
Grafic_afisare(Plot, "T_5_Gook_SCM_HF_raw");
dn_normat = Normare_Spectru(copy(dn), T_M);
Plot = Grafic_plot(E, [1.0 for i in eachindex(E)], "Fitul datelor experimentale ale spectrului neutronilor prompți cu o distribuție Maxwell", "", "E [MeV]", "Raportul spectrului neutronilor prompți la distribuția Maxwell", 0.0, 4.5, :log10, :identity, 10.0, :red);
Plot = Grafic_scatter(dn_normat.E, dn_normat.N, dn_normat.σN, :black, "Date experimentale ale lui Gook în SCM pentru HF", Plot);
Plot = Grafic_textbox(2.6, 4.1, latexstring("\$\\chi^2\$ = $(round(χ_sq, digits = 3)), \$\\mathrm{T_M}\$ = $(round(T_M, digits = 3)) MeV"), Plot);
Grafic_afisare(Plot, "T_5_Gook_SCM_HF_normat");