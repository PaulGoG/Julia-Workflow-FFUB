using Plots
using CSV
using DataFrames
using LaTeXStrings
using LsqFit

gr();

cd(@__DIR__) #Adauga path-ul la fisier

struct Separare
    Z
    A
    S
    σˢ
end

# Calculul energiei de separare a particulei cu A_part si Z_part din nucleu
function Energie_separare(librarie, A_part, Z_part)
    # Citim fisierul de tip CSV
    # Z|A|Simbol|D(KeV)|σᴰ(KeV) -> forma tabelului
    df = CSV.File(librarie; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame
    # Verifica daca perechea (A, Z) a particulei separate exista in libraria de date
    if isassigned(df.D[(df.A .== A_part) .& (df.Z .== Z_part)], 1)
        separare = Separare(Int[], Int[], Float64[], Float64[])
        # Citim defectul de masa pentru particula separata
        D_part = df.D[(df.A .== A_part) .& (df.Z .== Z_part)][1]
        σ_part = df.σ[(df.A .== A_part) .& (df.Z .== Z_part)][1]

        # Parcurgem toti radionuclizii din baza de date dupa A si dupa Z
        for i in (A_part + 1):maximum(df.A)
            Z_min = minimum(df.Z[df.A .== i])
            Z_max = maximum(df.Z[df.A .== i])
            for j in Z_min:Z_max 
                # Ne asiguram ca exista perechea (A-A_particula, Z-Z_particula)
                if j >= Z_part && isassigned(df.D[(df.A .== i - A_part) .& (df.Z .== j - Z_part)], 1)
                    D = df.D[(df.A .== i - A_part) .& (df.Z .== j - Z_part)][1]
                    σᴰ = df.σ[(df.A .== i - A_part) .& (df.Z .== j - Z_part)][1]

                    # Calculul valorilor
                    S = D + D_part - df.D[(df.A .== i) .& (df.Z .== j)][1]
                    σˢ = sqrt(σ_part^2 + σᴰ^2 + (df.σ[(df.A .== i) .& (df.Z .== j)][1])^2)
                    
                    # Scrierea valorilor in structura de vectori de vectori
                    # Valorile sunt calculate in MeVi!
                    push!(separare.S, S *1e-3)
                    push!(separare.σˢ, σˢ *1e-3)
                    push!(separare.A, i)
                    push!(separare.Z, j)
                end
            end
        end 
        return separare # Obiectul returnat de functie este un vector de vectori de tipul structului
    else 
        println("Perechea (A,Z) pentru particula pe care vrem sa o separam din nucleu nu exista in libraria de date $(librarie)") 
    end
end

# Stergem valorile negative ale energiei de separare acolo unde este cazul
function Valori_negative(separare)
    deleteat!(separare.A, findall(x -> x <0, separare.S))
    deleteat!(separare.Z, findall(x -> x <0, separare.S))
    deleteat!(separare.σˢ, findall(x -> x <0, separare.S))
    deleteat!(separare.S, findall(x -> x <0, separare.S))   
    return separare
end

# Reprezentari grafice
# Scatter simplu
function Grafic_simplu(separare, librarie, particula, scalare)
    upper_bound = maximum(separare.S)*scalare
    if minimum(separare.S) < 0 
        lower_bound = minimum(separare.S)*scalare
    else
        lower_bound = 0
    end
    plt = scatter(
    separare.A, 
    separare.S,  
    marker = :circle,
    markersize = 4, 
    markerstrokewidth = 0,
    ylims = (lower_bound, upper_bound),
    xlabel = L"\mathrm{A}", 
    ylabel = "Energia de separare a $(particula) [MeV]", 
    framestyle = :box,
    legend = :false,
    title = "Energia de separare a $(particula) în funcție de A, $(librarie[begin:end-4])",
    minorgrid = :true,
    mc = :blue, 
    size = (1000, 1000)
    )
    if lower_bound < 0
        hline!([0], ls = :dashdot, label = "")
    end
    display(plt)
    savefig(plt, "Grafice\\S_monocolor_$(librarie[begin:end-4]).png")
end
# Scatter impartit pe paritati
function Grafice_paritati_combinate(separare, librarie, particula, scalare)
    upper_bound = maximum(separare.S)*scalare
    if minimum(separare.S) < 0 
        lower_bound = minimum(separare.S)*scalare
    else
        lower_bound = 0
    end
    plt = scatter(
    separare.A[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 1)], 
    separare.S[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 1)],  
    marker = :xcross,
    markerstrokewidth = 0,
    markersize = 4, 
    ylims = (lower_bound, upper_bound),
    xlabel = L"\mathrm{A}", 
    ylabel = "Energia de separare a $(particula) [MeV]", 
    framestyle = :box,
    label = "par-par",
    title = "Energia de separare a $(particula) în funcție de A, $(librarie[begin:end-4])",
    minorgrid = :true,
    mc = :red, 
    size = (1000, 1000)
    )
    plt = scatter!(
    separare.A[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 0)], 
    separare.S[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 0)],  
    marker = :utriangle,
    markersize = 4, 
    markerstrokewidth = 0,
    xlabel = L"\mathrm{A}", 
    ylabel = "Energia de separare a $(particula) [MeV]", 
    framestyle = :box,
    label = "impar-impar",
    minorgrid = :true,
    mc = :green, 
    size = (1000, 1000)
    )
    plt = scatter!(
    separare.A[(iseven.(separare.A) .== 0)], 
    separare.S[(iseven.(separare.A) .== 0)],  
    marker = :star5,
    markersize = 4, 
    markerstrokewidth = 0,
    xlabel = L"\mathrm{A}", 
    ylabel = "Energia de separare a $(particula) [MeV]", 
    framestyle = :box,
    label = "A impar",
    minorgrid = :true,
    mc = :lightskyblue, 
    size = (1000, 1000)
    )     
    if lower_bound < 0
        hline!([0], ls = :dashdot, label = "")
    end
    display(plt)
    savefig(plt, "Grafice\\S_paritati_combinat_$(librarie[begin:end-4]).png")    
end
# Fit N/Z al Sₙ pentru toate elementele
function Grafic_fitare_simplu_neutron(librarie, scalare)
    separare = Valori_negative(Energie_separare(audi95, 1, 0))
    fitare(t, p) = p[1] .+ p[2].*t .+ p[3].*t.^2

    y = separare.S
    x = (separare.A - separare.Z) ./ (separare.Z)
    deleteat!(y, findall(z -> z>2.2, x))
    deleteat!(x, findall(z -> z>2.2, x))
    p0 = [1.0, 1.0, 1.0]
    fit = curve_fit(fitare, x, y, p0)
    A = round(fit.param[1], digits = 7)
    B = round(fit.param[2], digits = 7)
    C = round(fit.param[3], digits = 7)
    if sign(B) == 1
        semn1 = "+"
    elseif sign(B) == -1
        semn1 = "-"
    end
    if sign(C) == 1
        semn2 = "+"
    elseif sign(C) == -1
        semn2 = "-"
    end
    plt = scatter(
        x,
        y,
        marker = :circle,
        markersize = 4, 
        markerstrokewidth = 0,
        ylims = (0, maximum(y)*scalare),
        xlims = (minimum(x)*(2 - scalare), maximum(x)*scalare),
        xlabel = L"\frac{N}{Z}", 
        ylabel = "Energia de separare a neutronului [MeV]", 
        framestyle = :box,
        label = latexstring("\$S_{\\textrm{n}}\$ al tuturor nucleelor"),
        title = "Fitul energiei de separare a neutronului în funcție de N/Z, $(librarie[begin:end-4])",
        minorgrid = :true,
        mc = :blue, 
        size = (1000, 1000)
    )
    plt = plot!(
        sort(x),
        z -> A + B*z + C*z^2,
        color = :red,
        label = latexstring("\$ $(A) $(semn1) $(abs(B)) \\; \\frac{N}{Z} $(semn2) $(abs(C)) \\; \\left( \\frac{N}{Z} \\right)^2\$")
    )
    display(plt)
    savefig(plt, "Grafice\\S_fitare_simplu_neutron.png")
end
# Fit N/Z al Sₙ impartit pe layout pentru paritati
function Grafic_fitare_multiplu_neutron(librarie, scalare)
    separare = Valori_negative(Energie_separare(librarie, 1, 0))
    fitare(t, p) = p[1] .+ p[2].*t .+ p[3].*t.^2

    # p-p
    y = separare.S[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 1)]
    x = (separare.A[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 1)] .- separare.Z[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 1)]) ./ (separare.Z[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 1)])
    deleteat!(y, findall(z -> z>2.2, x))
    deleteat!(x, findall(z -> z>2.2, x))
    p0 = [1.0, 1.0, 1.0]
    fit = curve_fit(fitare, x, y, p0)
    A = round(fit.param[1], digits = 7)
    B = round(fit.param[2], digits = 7)
    C = round(fit.param[3], digits = 7)
    if sign(B) == 1
        semn1 = "+"
    elseif sign(B) == -1
        semn1 = "-"
    end
    if sign(C) == 1
        semn2 = "+"
    elseif sign(C) == -1
        semn2 = "-"
    end
    plt1 = scatter(
        x,
        y,
        marker = :circle,
        markersize = 4, 
        markerstrokewidth = 0,
        ylims = (0, maximum(y)*scalare),
        xlims = (minimum(x)*(2 - scalare), maximum(x)*scalare),
        xlabel = "", 
        ylabel = "Energia de separare a neutronului [MeV]", 
        framestyle = :box,
        label = "p-p",
        minorgrid = :true,
        mc = :lightblue, 
        size = (1000, 1000)
    )
    plt1 = plot!(
        sort(x),
        z -> A + B*z + C*z^2,
        color = :red,
        label = latexstring("\$ $(A) $(semn1) $(abs(B)) \\; \\frac{N}{Z} $(semn2) $(abs(C)) \\; \\left( \\frac{N}{Z} \\right)^2\$")
    )

    # i-i
    y = separare.S[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 0)]
    x = (separare.A[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 0)] .- separare.Z[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 0)]) ./ (separare.Z[(iseven.(separare.A) .== 1) .& (iseven.(separare.Z) .== 0)])
    deleteat!(y, findall(z -> z>2.2, x))
    deleteat!(x, findall(z -> z>2.2, x))
    p0 = [1.0, 1.0, 1.0]
    fit = curve_fit(fitare, x, y, p0)
    A = round(fit.param[1], digits = 7)
    B = round(fit.param[2], digits = 7)
    C = round(fit.param[3], digits = 7)
    if sign(B) == 1
        semn1 = "+"
    elseif sign(B) == -1
        semn1 = "-"
    end
    if sign(C) == 1
        semn2 = "+"
    elseif sign(C) == -1
        semn2 = "-"
    end
    plt2 = scatter(
        x,
        y,
        marker = :circle,
        markersize = 4, 
        markerstrokewidth = 0,
        ylims = (0, maximum(y)*scalare),
        xlims = (minimum(x)*(2 - scalare), maximum(x)*scalare),
        xlabel = "", 
        ylabel = "", 
        framestyle = :box,
        label = "i-i",
        minorgrid = :true,
        mc = :lightblue, 
        size = (1000, 1000)
    )
    plt2 = plot!(
        sort(x),
        z -> A + B*z + C*z^2,
        color = :red,
        label = latexstring("\$ $(A) $(semn1) $(abs(B)) \\; \\frac{N}{Z} $(semn2) $(abs(C)) \\; \\left( \\frac{N}{Z} \\right)^2\$")
    )

    # p-i
    y = separare.S[(iseven.(separare.A) .== 0) .& (iseven.(separare.Z) .== 1)]
    x = (separare.A[(iseven.(separare.A) .== 0) .& (iseven.(separare.Z) .== 1)] .- separare.Z[(iseven.(separare.A) .== 0) .& (iseven.(separare.Z) .== 1)]) ./ (separare.Z[(iseven.(separare.A) .== 0) .& (iseven.(separare.Z) .== 1)])
    deleteat!(y, findall(z -> z>2.2, x))
    deleteat!(x, findall(z -> z>2.2, x))
    p0 = [1.0, 1.0, 1.0]
    fit = curve_fit(fitare, x, y, p0)
    A = round(fit.param[1], digits = 7)
    B = round(fit.param[2], digits = 7)
    C = round(fit.param[3], digits = 7)
    if sign(B) == 1
        semn1 = "+"
    elseif sign(B) == -1
        semn1 = "-"
    end
    if sign(C) == 1
        semn2 = "+"
    elseif sign(C) == -1
        semn2 = "-"
    end
    plt3 = scatter(
        x,
        y,
        marker = :circle,
        markersize = 4, 
        markerstrokewidth = 0,
        ylims = (0, maximum(y)*scalare),
        xlims = (minimum(x)*(2 - scalare), maximum(x)*scalare),
        xlabel = L"\frac{N}{Z}", 
        ylabel = "Energia de separare a neutronului [MeV]", 
        framestyle = :box,
        label = "p-i",
        minorgrid = :true,
        mc = :lightblue, 
        size = (1000, 1000)
    )
    plt3 = plot!(
        sort(x),
        z -> A + B*z + C*z^2,
        color = :red,
        label = latexstring("\$ $(A) $(semn1) $(abs(B)) \\; \\frac{N}{Z} $(semn2) $(abs(C)) \\; \\left( \\frac{N}{Z} \\right)^2\$")
    )

    # i-p
    y = separare.S[(iseven.(separare.A) .== 0) .& (iseven.(separare.Z) .== 0)]
    x = (separare.A[(iseven.(separare.A) .== 0) .& (iseven.(separare.Z) .== 0)] .- separare.Z[(iseven.(separare.A) .== 0) .& (iseven.(separare.Z) .== 0)]) ./ (separare.Z[(iseven.(separare.A) .== 0) .& (iseven.(separare.Z) .== 0)])
    deleteat!(y, findall(z -> z>2.2, x))
    deleteat!(x, findall(z -> z>2.2, x))
    p0 = [1.0, 1.0, 1.0]
    fit = curve_fit(fitare, x, y, p0)
    A = round(fit.param[1], digits = 7)
    B = round(fit.param[2], digits = 7)
    C = round(fit.param[3], digits = 7)
    if sign(B) == 1
        semn1 = "+"
    elseif sign(B) == -1
        semn1 = "-"
    end
    if sign(C) == 1
        semn2 = "+"
    elseif sign(C) == -1
        semn2 = "-"
    end
    plt4 = scatter(
        x,
        y,
        marker = :circle,
        markersize = 4, 
        markerstrokewidth = 0,
        ylims = (0, maximum(y)*scalare),
        xlims = (minimum(x)*(2 - scalare), maximum(x)*scalare),
        xlabel = L"\frac{N}{Z}", 
        ylabel = "", 
        framestyle = :box,
        label = "i-p",
        minorgrid = :true,
        mc = :lightblue, 
        size = (1000, 1000)
    )
    plt4 = plot!(
        sort(x),
        z -> A + B*z + C*z^2,
        color = :red,
        label = latexstring("\$ $(A) $(semn1) $(abs(B)) \\; \\frac{N}{Z} $(semn2) $(abs(C)) \\; \\left( \\frac{N}{Z} \\right)^2\$")
    )
    plt = plot(plt1, plt2, plt3, plt4, layout = (2, 2))
    display(plt)
    savefig(plt, "Grafice\\S_fitare_multiplu_neutron.png")
end
# Grafice S in functie de N si Z
function Grafic_dublu_NZ(separare, librarie, particula, scalare)
    upper_bound = maximum(separare.S)*scalare
    if minimum(separare.S) < 0 
        lower_bound = minimum(separare.S)*scalare
    else
        lower_bound = 0
    end
    plt1 = scatter(
    separare.Z, 
    separare.S,  
    marker = :circle,
    markersize = 4, 
    markerstrokewidth = 0,
    ylims = (lower_bound, upper_bound),
    xlabel = L"\mathrm{Z}", 
    ylabel = "Energia de separare a $(particula) [MeV]", 
    framestyle = :box,
    legend = :false,
    title = "Energia de separare a $(particula) în funcție de Z și N, $(librarie[begin:end-4])",
    minorgrid = :true,
    mc = :purple, 
    size = (1000, 1000)
    )
    if lower_bound < 0
        hline!([0], ls = :dashdot, label = "")
    end
    plt2 = scatter(
    (separare.A .- separare.Z), 
    separare.S,  
    marker = :circle,
    markersize = 4, 
    markerstrokewidth = 0,
    ylims = (lower_bound, upper_bound),
    xlabel = L"\mathrm{N}", 
    ylabel = "Energia de separare a $(particula) [MeV]", 
    framestyle = :box,
    legend = :false,
    minorgrid = :true,
    mc = :green, 
    size = (1000, 1000)
    )
    if lower_bound < 0
        hline!([0], ls = :dashdot, label = "")
    end
    plt = plot(plt1, plt2, layout = (2,1))
    display(plt)
    savefig(plt, "Grafice\\S_dublu_NZ_$(librarie[begin:end-4]).png")
end

# Apelarea functiilor definite
audi95 = "AUDI95.csv"
audi21 = "AUDI2021.csv"
scalare = 1.05

separare_n_95 = Valori_negative(Energie_separare(audi95, 1, 0))
separare_p_95 = Valori_negative(Energie_separare(audi95, 1, 1))
separare_α_95 = Energie_separare(audi95, 4, 2)
separare_d_95 = Valori_negative(Energie_separare(audi95, 2, 1))

Grafice_paritati_combinate(separare_n_95, audi95, L"^1_0n", scalare)
Grafice_paritati_combinate(separare_p_95, audi95, L"^1_1p", scalare)
Grafic_simplu(separare_α_95, audi95, L"^4_2\alpha", scalare)
Grafic_simplu(separare_d_95, audi95, L"^2_1H", scalare)
Grafic_fitare_simplu_neutron(audi95, scalare)
Grafic_fitare_multiplu_neutron(audi95, scalare)
Grafic_dublu_NZ(separare_n_95, audi95, L"^1_0n", scalare)
Grafic_dublu_NZ(separare_p_95, audi95, L"^1_1p", scalare)