using Plots
using CSV
using DataFrames
using LaTeXStrings

gr();

cd(@__DIR__) #Adauga path-ul la fisier

struct Radionuclid
    Z
    A
    B
    σᴮ
end

struct Diferente
    Z
    A
    εᴰ
    σᴰ
end

function Energie_medie_legatura(librarie)
    radionuclid = Radionuclid(Int[], Int[], Float64[], Float64[])
    # Citim fisierul de tip CSV
    # Z|A|Simbol|D(KeV)|σᴰ(KeV) -> forma tabelului
    df = CSV.File(librarie; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame

    # Citim defectele de masa si incertitudinile pentru neutron si proton
    Dᵖ = df.D[(df.A .== 1) .& (df.Z .== 1)][1]
    σᵖ = df.σ[(df.A .== 1) .& (df.Z .== 1)][1]
    Dⁿ = df.D[(df.A .== 1) .& (df.Z .== 0)][1]
    σⁿ = df.σ[(df.A .== 1) .& (df.Z .== 0)][1]

    # Parcurgem toti radionuclizii din baza de date dupa A si dupa Z
    for i in 2:maximum(df.A)
        Z_min = minimum(df.Z[df.A .== i])
        Z_max = maximum(df.Z[df.A .== i])
        for j in Z_min:Z_max 
            # Selectare rubrici de interes pentru elementul (A,Z)
            D = df.D[(df.A .== i) .& (df.Z .== j)][1]
            σᴰ = df.σ[(df.A .== i) .& (df.Z .== j)][1]

            # Calculul valorilor
            W = j*Dᵖ + (i-j)*Dⁿ - D
            σᵂ = sqrt(j^2 * σᵖ^2 + (i-j)^2 * σⁿ^2 + σᴰ^2)
            
            # Scrierea valorilor in structura de vectori de vectori
            # Valorile energiilor sunt calculate in MeVi
            push!(radionuclid.B, W/i *1e-3)
            push!(radionuclid.σᴮ, σᵂ/i *1e-3)
            push!(radionuclid.A, i)
            push!(radionuclid.Z, j)
        end
    end 

    # Figuram doar valorile mai mari de 3 MeV pentru scalarea buna a axelor
    deleteat!(radionuclid.A, findall(x -> x <=3, radionuclid.B))
    deleteat!(radionuclid.Z, findall(x -> x <=3, radionuclid.B))
    deleteat!(radionuclid.σᴮ, findall(x -> x <=3, radionuclid.B))
    deleteat!(radionuclid.B, findall(x -> x <= 3, radionuclid.B))

    return radionuclid # Obiectul returnat de functie este un vector de vectori de tipul structului
end

function Diferente_relative(librarie1, librarie2, trunchiere_sup, trunchiere_inf)
    diferente = Diferente(Int[], Int[], Float64[], Float64[])
    # Citim fisierele
    # Z|A|Simbol|D(KeV)|σᴰ(KeV) -> forma tabelului
    df1 = CSV.File(librarie1; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame
    df2 = CSV.File(librarie2; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame

    # Indexam parcurgerea librariilor de date dupa prima, si impunem existenta campului de date in cea de-a doua
    for i in 1:maximum(df1.A)
            Z_min = minimum(df1.Z[df1.A .== i])
            Z_max = maximum(df1.Z[df1.A .== i])
            for j in Z_min:Z_max 
                # Verifica daca elementul (A,Z) exista in a doua librarie de date
                if isassigned(df2.Z[df2.A .== i], j - Z_min + 1) 
                    D1 = abs(df1.D[(df1.A .== i) .& (df1.Z .== j)][1])
                    σᴰ1 = df1.σ[(df1.A .== i) .& (df1.Z .== j)][1]
                    D2 = abs(df2.D[(df2.A .== i) .& (df2.Z .== j)][1])
                    σᴰ2 = df2.σ[(df2.A .== i) .& (df2.Z .== j)][1]

                    if D1 !=0 && D2 !=0                         
                        εᴰ = 100 * abs(D1-D2)/D1 # Diferenta relativa data in procente
                        if εᴰ <= trunchiere_sup && εᴰ >= trunchiere_inf # Ignoram valorile peste un procent dat pentru buna scalare a graficului
                            σᴰ = 100 * sqrt((1 + D2/D1^2)^2 * σᴰ1^2 + (1 - 1/D1)^2 * σᴰ2^2) # Incertitudinea
                                
                            push!(diferente.εᴰ, εᴰ)
                            push!(diferente.σᴰ, σᴰ)
                            push!(diferente.A, i)
                            push!(diferente.Z, j)
                        end
                    end
                end
            end
    end 
    return diferente
end

# Reprezentari grafice
# Scatter simplu
function Grafic_simplu(radionuclid, librarie)
    plt = scatter(
    radionuclid.A, 
    radionuclid.B,  
    marker = :circle,
    markersize = 3, 
    markerstrokewidth = 0,
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z) = \\frac{W(A,Z)}{A}\$  [MeV]"), 
    framestyle = :box,
    legend = :false,
    title = "Energia medie de legătură per nucleon, $(librarie[begin:end-4])",
    minorgrid = :true,
    mc = :limegreen, 
    size = (1000, 1000)
    )
    display(plt)
    savefig(plt, "Grafice\\B_toti_radionuclizii_$(librarie[begin:end-4]).png")
end
# Scatter cu error bars
function Grafic_errorb(radionuclid, librarie)
    plt1 = scatter(
    radionuclid.A, 
    radionuclid.B, 
    yerr = radionuclid.σᴮ, 
    marker = :xcross,
    markersize = 3, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z) = \\frac{W(A,Z)}{A}\$  [MeV]"), 
    framestyle = :box,
    legend = :false,
    title = "Energia medie de legătură per nucleon, $(librarie[begin:end-4])",
    minorgrid = :true,
    msc = :red,
    mc = :aqua, 
    size = (1000, 1000)
    )
    annotate!(150, 6.5, L"\sigma_{B(A,Z)} = \frac{1}{A} \sqrt{Z^2 \sigma_{D_p}^2 + (A-Z)^2 \sigma_{D_n}^2 + \sigma_{D(A,Z)}^2}")
    display(plt1)
    savefig(plt1, "Grafice\\B_toti_radionuclizii_errb_$(librarie[begin:end-4]).png")
end
# Scatter comparativ intre 2 biblioteci de date
function Grafic_comparativ_lib(radionuclid_1, radionuclid_2, librarie_1, librarie_2)
    plt1 = scatter(
    radionuclid_1.A, 
    radionuclid_1.B,  
    marker = :circle,
    markersize = 3,
    markerstrokewidth = 0, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z)\$ [MeV]"), 
    framestyle = :box,
    legend = :false,
    title = "Energia medie de legătură per nucleon, $(librarie_1[begin:end-4])",
    minorgrid = :true,
    mc = :red, 
    size = (900, 900)
    )
    plt2 = scatter(
    radionuclid_2.A, 
    radionuclid_2.B, 
    marker = :circle,
    markersize = 3,
    markerstrokewidth = 0, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z)\$ [MeV]"), 
    framestyle = :box,
    legend = :false,
    title = "Energia medie de legătură per nucleon, $(librarie_2[begin:end-4])",
    minorgrid = :true,
    mc = :green, 
    size = (900, 900)
    )
    plt = plot(plt1, plt2, layout = (2, 1))
    display(plt)
    savefig(plt, "Grafice\\B_comparat_$(librarie_1[begin:end-4])_$(librarie_2[begin:end-4]).png")
end
# Scatter in functie de paritatea A si Z
function Grafice_paritati_combinate(radionuclid, librarie)
    plt = scatter(
    radionuclid.A[(iseven.(radionuclid.A) .== 1) .& (iseven.(radionuclid.Z) .== 1)], 
    radionuclid.B[(iseven.(radionuclid.A) .== 1) .& (iseven.(radionuclid.Z) .== 1)],  
    marker = :xcross,
    markerstrokewidth = 0,
    markersize = 3, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z)\$  [MeV]"), 
    framestyle = :box,
    label = "par-par",
    title = "Energia medie de legătură per nucleon, $(librarie[begin:end-4])",
    minorgrid = :true,
    mc = :red, 
    size = (1000, 1000)
    )
    plt = scatter!(
    radionuclid.A[(iseven.(radionuclid.A) .== 1) .& (iseven.(radionuclid.Z) .== 0)], 
    radionuclid.B[(iseven.(radionuclid.A) .== 1) .& (iseven.(radionuclid.Z) .== 0)],  
    marker = :utriangle,
    markersize = 3, 
    markerstrokewidth = 0,
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z)\$  [MeV]"), 
    framestyle = :box,
    label = "impar-impar",
    minorgrid = :true,
    mc = :green, 
    size = (1000, 1000)
    )
    plt = scatter!(
    radionuclid.A[(iseven.(radionuclid.A) .== 0)], 
    radionuclid.B[(iseven.(radionuclid.A) .== 0)],  
    marker = :star5,
    markersize = 3, 
    markerstrokewidth = 0,
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z)\$  [MeV]"), 
    framestyle = :box,
    label = "A impar",
    minorgrid = :true,
    mc = :lightskyblue, 
    size = (1000, 1000)
    )     
    display(plt)
    savefig(plt, "Grafice\\B_paritati_combinat_$(librarie[begin:end-4]).png")    
end
# Scatter cu layout 3 linii, o coloana
function Grafice_paritati_comparativ(radionuclid, librarie)
    plt1 = scatter(
    radionuclid.A[(iseven.(radionuclid.A) .== 1) .& (iseven.(radionuclid.Z) .== 1)], 
    radionuclid.B[(iseven.(radionuclid.A) .== 1) .& (iseven.(radionuclid.Z) .== 1)],  
    marker = :circle,
    markerstrokewidth = 0,
    markersize = 3, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z)\$  [MeV]"), 
    framestyle = :box,
    label = "par-par",
    title = "Energia medie de legătură per nucleon, $(librarie[begin:end-4])",
    minorgrid = :true,
    mc = :red, 
    size = (900, 900)
    )
    plt2 = scatter(
    radionuclid.A[(iseven.(radionuclid.A) .== 1) .& (iseven.(radionuclid.Z) .== 0)], 
    radionuclid.B[(iseven.(radionuclid.A) .== 1) .& (iseven.(radionuclid.Z) .== 0)],  
    marker = :utriangle,
    markersize = 3, 
    markerstrokewidth = 0,
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z)\$  [MeV]"), 
    framestyle = :box,
    label = "impar-impar",
    minorgrid = :true,
    mc = :green, 
    size = (900, 900)
    )
    plt3 = scatter(
    radionuclid.A[(iseven.(radionuclid.A) .== 0)], 
    radionuclid.B[(iseven.(radionuclid.A) .== 0)],  
    marker = :star5,
    markersize = 3, 
    markerstrokewidth = 0,
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z)\$  [MeV]"), 
    framestyle = :box,
    label = "A impar",
    minorgrid = :true,
    mc = :lightskyblue, 
    size = (900, 900)
    )   
    plt = plot(plt1, plt2, plt3, layout = (3, 1))  
    display(plt)
    savefig(plt, "Grafice\\B_paritati_comparativ_$(librarie[begin:end-4]).png")    
end
# Scatter diferente relative intre masele a doua librarii de date
function Grafic_diferente(diferente, librarie1, librarie2, trunchiere_sup, trunchiere_inf)
    plt1 = scatter(
    diferente.A, 
    diferente.εᴰ, 
    #yerr = diferente.σᴰ, 
    marker = :xcross,
    markersize = 3, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$\\varepsilon_{\\textrm{D}} = \\frac{|D_{\\textrm{$(librarie1[begin:end-4])}} - D_{\\textrm{$(librarie2[begin:end-4])}}|}{D_{\\textrm{$(librarie1[begin:end-4])}}}\$ exprimat în procente"), 
    framestyle = :box,
    legend = :false,
    title = "Diferentele relative ale defectelor de masă între librăriile $(librarie1[begin:end-4]) și $(librarie2[begin:end-4])",
    minorgrid = :true,
    msc = :red,
    mc = :blue, 
    size = (1000, 1000)
    )
    annotate!(diferente.A[Int(ceil(length(diferente.A)/2))], Int(ceil(trunchiere_sup*0.8)), latexstring("\$\\varepsilon_{\\textrm{D}} \\in [$(trunchiere_inf)\\%, $(trunchiere_sup)\\%]\$"))
    #annotate!(diferente.A[Int(ceil(length(diferente.A)/2))], Int(ceil(trunchiere_sup*0.6)), latexstring("\$\\sigma_{\\varepsilon_{\\textrm{D}}} = 100 \\sqrt{(1 + \\frac{D_{\\textrm{$(librarie2[begin:end-4])}}}{D^2_{\\textrm{$(librarie1[begin:end-4])}}})^2 \\sigma^2_{D_{\\textrm{$(librarie1[begin:end-4])}}} + (1 - \\frac{1}{D_{\\textrm{$(librarie1[begin:end-4])}}})^2 \\sigma^2_{D_{\\textrm{$(librarie2[begin:end-4])}}}}\$"))
    display(plt1)
    savefig(plt1, "Grafice\\Diferente_relative_defecte_mase_errb_$(librarie1[begin:end-4])_$(librarie2[begin:end-4])_trunch_$(trunchiere_inf)_$(trunchiere_sup).png")
end

# Apelarea functiilor definite
audi95 = "AUDI95.csv"
audi21 = "AUDI2021.csv"
trunchiere_sup = 10
trunchiere_inf = 0.05
radionuclid_95 = Energie_medie_legatura(audi95)
radionuclid_21 = Energie_medie_legatura(audi21)
diferente = Diferente_relative(audi95, audi21, trunchiere_sup, trunchiere_inf)

Grafic_simplu(radionuclid_95, audi95)
Grafic_simplu(radionuclid_21, audi21)
Grafic_errorb(radionuclid_95, audi95)
Grafic_errorb(radionuclid_21, audi21)
Grafic_comparativ_lib(radionuclid_95, radionuclid_21, audi95, audi21)
Grafice_paritati_combinate(radionuclid_95, audi95)
Grafice_paritati_combinate(radionuclid_21, audi21)
Grafice_paritati_comparativ(radionuclid_95, audi95)
Grafice_paritati_comparativ(radionuclid_21, audi21)
Grafic_diferente(diferente, audi95, audi21, trunchiere_sup, trunchiere_inf)
Grafic_diferente(diferente, audi21, audi95, trunchiere_sup, trunchiere_inf)