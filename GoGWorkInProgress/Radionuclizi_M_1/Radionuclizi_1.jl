using Plots
using CSV
using DataFrames
using LaTeXStrings

#plotlyjs();

cd(@__DIR__) #Adauga path-ul la fisier

struct Radionuclid
    Z
    A
    B
    σᴮ
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

    for i in 2:maximum(df.A)
        Z_min = minimum(df.Z[df.A .== i])
        Z_max = maximum(df.Z[df.A .== i])
        for j in Z_min:Z_max # Parcurgem toti radionuclizii din baza de date dupa A si dupa Z
            D = df.D[(df.A .== i) .& (df.Z .== j)][1]
            σᴰ = df.σ[(df.A .== i) .& (df.Z .== j)][1]

            W = j*Dᵖ + (i-j)*Dⁿ - D
            σᵂ = sqrt(j^2 * σᵖ^2 + (i-j)^2 * σⁿ^2 + σᴰ^2)

            push!(radionuclid.B, W/i *1e-3)
            push!(radionuclid.σᴮ, σᵂ/i *1e-3)
            push!(radionuclid.A, i)
            push!(radionuclid.Z, j)
        end
    end 

    # Figuram doar valorile mai mari de 3 MeV pentru scalare buna a axelor
    deleteat!(radionuclid.A, findall(x -> x <=3, radionuclid.B))
    deleteat!(radionuclid.σᴮ, findall(x -> x <=3, radionuclid.B))
    deleteat!(radionuclid.B, findall(x -> x <= 3, radionuclid.B))

    # Reprezentari grafice
    plt1 = scatter(
    radionuclid.A, 
    radionuclid.B,  
    marker = :xcross,
    markersize = 3, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z) = \\frac{W(A,Z)}{A}\$  [MeV]"), 
    framestyle = :box,
    legend = :false,
    title = "Energia medie de legătură per nucleon",
    minorgrid = :true,
    mc = :blue, 
    size = (1000, 1000)
    )
    display(plt1)
    savefig(plt1, "Grafice\\B_toti_radionuclizii_$(librarie[begin:end-4]).png")

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
    title = "Energia medie de legătură per nucleon",
    minorgrid = :true,
    msc = :olivedrab,
    mc = :cyan, 
    size = (1000, 1000)
    )
    annotate!(150, 6.5, L"\sigma_{B(A,Z)} = \frac{1}{A} \sqrt{Z^2 \sigma_{D_p}^2 + (A-Z)^2 \sigma_{D_n}^2 + \sigma_{D(A,Z)}^2}")
    savefig(plt1, "Grafice\\B_toti_radionuclizii_errb_$(librarie[begin:end-4]).png")

    return radionuclid # Obiectul returnat de functie este un vector de vectori 
end

radionuclid = Energie_medie_legatura("AUDI95.csv")

    plt1 = scatter(
    radionuclid.A, 
    radionuclid.B,  
    marker = :xcross,
    markersize = 3, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z) = \\frac{W(A,Z)}{A}\$  [MeV]"), 
    framestyle = :box,
    legend = :false,
    title = "Energia medie de legătură per nucleon",
    minorgrid = :true,
    mc = :blue, 
    size = (900, 900)
    )

    plt2 = scatter(
    radionuclid.A, 
    radionuclid.B, 
    yerr = radionuclid.σᴮ, 
    marker = :xcross,
    markersize = 3, 
    xlabel = L"\mathrm{A}", 
    ylabel = latexstring("\$B(A,Z) = \\frac{W(A,Z)}{A}\$  [MeV]"), 
    framestyle = :box,
    legend = :false,
    title = "Energia medie de legătură per nucleon",
    minorgrid = :true,
    msc = :olivedrab,
    mc = :cyan, 
    size = (900, 900)
    )
    annotate!(150, 6.5, L"\sigma_{B(A,Z)} = \frac{1}{A} \sqrt{Z^2 \sigma_{D_p}^2 + (A-Z)^2 \sigma_{D_n}^2 + \sigma_{D(A,Z)}^2}")
    plt = plot(plt1, plt2, layout = (2, 1))
    savefig(plt1, "Grafice\\B_toti_radionuclizii_errb_$(librarie[begin:end-4]).png")