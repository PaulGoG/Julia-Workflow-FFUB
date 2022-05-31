using Plots
using CSV
using DataFrames
using LaTeXStrings

# Calculul corectiei de paturi si realizarea graficelor aferente temei 5

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Definirea obiectului in care vom stoca datele calculate
struct Paturi
    Z
    A
    W_exp
    W_LDM
    δW_0
    δW
    P_a
end

# Citirea datelor global din libraria de date
librarie = "AUDI95.csv"
df = CSV.File(librarie; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame
dm = CSV.File("MOLLER.csv"; delim=' ', ignorerepeated=true, header=["Z", "A", "δW"]) |> DataFrame

# Termenii parametrizarii Pearson in MeV
a_v = 15.65;
a_s = 17.63;
a_c = 0.864/1.233;
a_s1 = 27.72;
a_s2 = 25.60;

# Citirea defectelor de masa pentru proton si neutron in MeV
Dᵖ = df.D[(df.A .== 1) .& (df.Z .== 1)][1]*1e-3;
Dⁿ = df.D[(df.A .== 1) .& (df.Z .== 0)][1]*1e-3;

# Calculul punctual al W_exp folosind excesele de masa
function Calcul_W_exp(A, Z)
    D = df.D[(df.A .== A) .& (df.Z .== Z)][1]*1e-3
    return Z*Dᵖ + (A-Z)*Dⁿ - D
end

# Calculul punctual al W_LDM fara termenul de imperechere conform parametrizarii Pearson
function Calcul_W_LDM(A, Z)
    a_sim = A*(a_s1 - a_s2* A^(-1/3))
    η = (A - 2*Z)/A
    return a_v*A - a_s*A^(2/3) - a_c*Z^2 *A^(-1/3) - a_sim*η^2
end

# Calculul punctual al termenului de imperechere, existenta D-urilor NU este verificata aici!
function Calcul_P_a(A, Z)
    D_plus = df.D[(df.A .== A+2) .& (df.Z .== Z+1)][1]*1e-3
    D_minus = df.D[(df.A .== A-2) .& (df.Z .== Z-1)][1]*1e-3
    D = df.D[(df.A .== A) .& (df.Z .== Z)][1]*1e-3
    return 0.5*(D_plus - 2*D + D_minus)
end

# Calculul punctual al corectiei de paturi fara termenul de imperechere
function Calcul_δW_0(A, Z)
    return Calcul_W_LDM(A, Z) - Calcul_W_exp(A, Z)
end

# Calculul punctual al corectiei de paturi finale
function Calcul_δW(A, Z)
    return Calcul_δW_0(A, Z) + 0.5*Calcul_P_a(A, Z)
end

# Functia finala de constructie a obiectului cu ajutorul caruia vom realiza reprezentarile grafice
function Constructie_struct()
    paturi = Paturi(Int[], Int[], Float64[], Float64[], Float64[], Float64[], Float64[])
    # Parcurgem toti radionuclizii din baza de date dupa A >= 3 si dupa Z
    for i in 3:maximum(df.A)
        Z_min = minimum(df.Z[df.A .== i])
        Z_max = maximum(df.Z[df.A .== i])
        for j in Z_min:Z_max
            # Impunere conditii de existenta ale elementelor din calculul P_a
            if isassigned(df.D[(df.A .== i-2) .& (df.Z .== j-1)], 1) && isassigned(df.D[(df.A .== i+2) .& (df.Z .== j+1)], 1)
                push!(paturi.A, i)
                push!(paturi.Z, j)
                push!(paturi.W_exp, Calcul_W_exp(i, j))
                push!(paturi.W_LDM, Calcul_W_LDM(i, j))
                push!(paturi.P_a, Calcul_P_a(i, j))
                push!(paturi.δW_0, Calcul_δW_0(i, j))
                push!(paturi.δW, Calcul_δW(i, j))
            end
        end
    end 
    return paturi
end

# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
# Comparare corectii de paturi cu si fara termen de imperechere
function Grafic_comparare_corectii(paturi)
    plt = scatter(
        paturi.A, 
        paturi.δW_0,  
        marker = :circle,
        markersize = 4, 
        markerstrokewidth = 0,
        xlabel = L"\mathrm{A}", 
        ylabel = latexstring("Corecțiile de pături \$\\delta W\$  [MeV]"), 
        framestyle = :box,
        label = latexstring("Fără termenul de împerechere \$\\mathrm{P_a}\$"),
        title = "Corecțiile de pături bazate pe parametrizarea Pearson a modelului picătura de lichid",
        minorgrid = :false,
        ylims = (-15, 10),
        yticks = [-15, -10, -5, 0, 5, 10],
        xlims = (0.0, 280.0),
        mc = :red, 
        size = (1000, 1000)
    )
    plt = scatter!(
        paturi.A, 
        paturi.δW,
        marker = :utriangle,
        markersize = 4, 
        markerstrokewidth = 0,
        mc = :lime,
        label = latexstring("Cu termenul de împerechere \$\\mathrm{P_a}\$"),
        size = (1000, 1000)
    )
    plt = hline!([0], label = "")
    display(plt)
    #savefig(plt, "Grafice\\Corectii_comparate.png") 
    #savefig(plt, "Grafice/Corectii_comparate.png")   
end
# Comparare corectie de paturi calculata vs Moller si Nix
function Grafic_Moller_Nix(paturi)
    plt = scatter(
        paturi.A, 
        paturi.δW,  
        marker = :circle,
        markersize = 4, 
        markerstrokewidth = 0,
        xlabel = L"\mathrm{A}", 
        ylabel = latexstring("Corecțiile de pături \$\\delta W\$  [MeV]"), 
        framestyle = :box,
        label = "Calculate cu parametrizarea Pearson a modelului picătura de lichid",
        title = "Comparație între rezultatele calculate cu \$\\mathrm{W_{LDM}}\$ și cele obținute de Moller și Nix",
        minorgrid = :false,
        xlims = (0.0, 280.0),
        ylims = (-20, 10),
        yticks = [-20, -15, -10, -5, 0, 5, 10],
        mc = :red, 
        size = (1000, 1000)
    )
    # Alegem din Moller si Nix doar radionuclizii care au corespondent calculat de noi
    y = zeros(0)
    x = zeros(0)
    for i in minimum(dm.A):maximum(dm.A)
        Z_min = minimum(dm.Z[dm.A .== i])
        Z_max = maximum(dm.Z[dm.A .== i])
        for j in Z_min:Z_max
            if isassigned(paturi.A[(paturi.A .== i) .& (paturi.Z .== j)], 1)
                push!(y, dm.δW[(dm.A .== i) .& (dm.Z .== j)][1])
                push!(x, dm.A[(dm.A .== i) .& (dm.Z .== j)][1])
            end
        end
    end
    plt = scatter!(
        x, 
        y,
        marker = :utriangle,
        markersize = 4, 
        markerstrokewidth = 0,
        mc = :blue,
        label = "Obținute de Moller și Nix (doar cele care au corespondent calculat)",
        size = (1000, 1000)
    )
    plt = hline!([0], label = "")
    display(plt)
    #savefig(plt, "Grafice\\Corectii_comparate_Moller_Nix.png")  
    #savefig(plt, "Grafice/Corectii_comparate_Moller_Nix.png")
end
# Reprezentari grafice termen perechi
function Grafice_perechi_paritati(paturi)
    # p-p
    y = paturi.P_a[(iseven.(paturi.A) .== 1) .& (iseven.(paturi.Z) .== 1)]
    x = paturi.A[(iseven.(paturi.A) .== 1) .& (iseven.(paturi.Z) .== 1)]
    plt = scatter(
        x,
        y,
        marker = :circle,
        markersize = 4, 
        markerstrokewidth = 0,
        title = "Termenul de împerechere Pearson defalcat pe parități",
        xlabel = L"\mathrm{A}", 
        ylabel = latexstring("\$\\mathrm{P_a}\$ împerechere Pearson [MeV]"), 
        ylims = (-8.5, 8.5),
        xlims = (0.0, 280.0),
        yticks = [-8, -6, -4, -2, 0, 2, 4, 6, 8],
        framestyle = :box,
        label = "par-par",
        minorgrid = :false,
        mc = :red, 
        size = (1000, 1000)
    )
    
    # i-i
    y = paturi.P_a[(iseven.(paturi.A) .== 1) .& (iseven.(paturi.Z) .== 0)]
    x = paturi.A[(iseven.(paturi.A) .== 1) .& (iseven.(paturi.Z) .== 0)]
    plt = scatter!(
        x,
        y,
        marker = :utriangle,
        markersize = 4, 
        markerstrokewidth = 0,
        label = "impar-impar",
        mc = :green, 
        size = (1000, 1000)
    )
    
    # p-i
    y = paturi.P_a[(iseven.(paturi.A) .== 0) .& (iseven.(paturi.Z) .== 1)]
    x = paturi.A[(iseven.(paturi.A) .== 0) .& (iseven.(paturi.Z) .== 1)]
    plt = scatter!(
        x,
        y,
        marker = :xcross,
        markersize = 4, 
        markerstrokewidth = 0,
        label = "par-impar",
        minorgrid = :true,
        mc = :lightblue, 
        size = (1000, 1000)
    )
    
    # i-p
    y = paturi.P_a[(iseven.(paturi.A) .== 0) .& (iseven.(paturi.Z) .== 0)]
    x = paturi.A[(iseven.(paturi.A) .== 0) .& (iseven.(paturi.Z) .== 0)]
    plt = scatter!(
        x,
        y,
        marker = :star5,
        markersize = 4, 
        markerstrokewidth = 0,
        label = "impar-par",
        minorgrid = :true,
        mc = :purple, 
        size = (1000, 1000)
    )
    plt = hline!([0], label = "")    
    display(plt)
    #savefig(plt, "Grafice\\Termen_perechi_Pearson.png")
    #savefig(plt, "Grafice/Termen_perechi_Pearson.png")
end

# Reprezentare grafica B = W/A 
function Grafic_B(paturi)
    y = paturi.W_exp./paturi.A
    plt = scatter(
        paturi.A, 
        y,  
        marker = :circle,
        markersize = 3, 
        markerstrokewidth = 0,
        xlabel = L"\mathrm{A}", 
        ylabel = "Energia medie de legătură per nucleon B(A, Z) [MeV]", 
        framestyle = :box,
        label = L"\mathrm{W}_{\mathrm{exp}}",
        title = latexstring("Comparație între valorile B(A, Z) calculate cu \$\\mathrm{W_{exp}}\$ și cele calculate cu \$\\mathrm{W_{LDM}}\$"),
        minorgrid = :true,
        xlims = (0.0, 280.0),
        ylims = (3.0, 9.0),
        mc = :lime, 
        size = (1000, 1000)
    )
    y = paturi.W_LDM./paturi.A
    plt = scatter!(
        paturi.A, 
        y,  
        marker = :xcross,
        markersize = 3, 
        markerstrokewidth = 0,
        label = L"\mathrm{W}_{\mathrm{LDM}}",
        mc = :lightblue, 
        size = (1000, 1000)
    ) 
    display(plt)
    #savefig(plt, "Grafice\\B_comparat_modele.png")  
    #savefig(plt, "Grafice/B_comparat_modele.png") 
end

# Apelarea functiilor definite pentru executia programului
paturi = Constructie_struct();

Grafice_perechi_paritati(paturi)
Grafic_comparare_corectii(paturi)
Grafic_Moller_Nix(paturi)
Grafic_B(paturi)