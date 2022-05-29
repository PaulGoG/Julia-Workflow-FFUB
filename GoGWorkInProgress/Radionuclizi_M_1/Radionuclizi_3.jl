using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul al energiei de imperechere & reprezentarile grafice asociate temei 3

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Definim obiectul in care stocam datele calculate
struct Pairing
    Z
    A
    P
    σᴾ
end

librarie = "AUDI95.csv"
df = CSV.File(librarie; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame

# Calculul energiei de separare a particulei cu A_part si Z_part din nucleul cu A si Z
function Energie_separare(A_part, Z_part, A, Z)
    # Verificarea existentei radionuclizilor folositi in libraria de date
    if isassigned(df.D[(df.A .== A_part) .& (df.Z .== Z_part)], 1) && isassigned(df.D[(df.A .== A) .& (df.Z .== Z)], 1) && isassigned(df.D[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)], 1)
        D_part = df.D[(df.A .== A_part) .& (df.Z .== Z_part)][1]
        #σ_part = df.σ[(df.A .== A_part) .& (df.Z .== Z_part)][1]
    
        D = df.D[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)][1]
        #σᴰ = df.σ[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)][1]    
    
        S = (D + D_part - df.D[(df.A .== A) .& (df.Z .== Z)][1])*1e-3
        #σˢ = sqrt(σ_part^2 + σᴰ^2 + (df.σ[(df.A .== A) .& (df.Z .== Z)][1])^2)*1e-3

        # Obiectul returnat este un vector cu variabila booleana pe prima pozitie si valoarea lui S pe a doua
        return [true, S]
    else 
        return [false, 0]
    end 
end

# Calculul energiei de imperechere pentru toti radionuclizii din libraria de date
function Energie_imperechere(formula)
    pairing = Pairing(Int[], Int[], Float64[], Float64[])

    for i in 1:maximum(df.A)
        Z_min = minimum(df.Z[df.A .== i])
        Z_max = maximum(df.Z[df.A .== i])
        for j in Z_min:Z_max
            Sn1 = Energie_separare(1,0,i+1,j)
            Sn2 = Energie_separare(1,0,i,j)
            Sn3 = Energie_separare(1,0,i-1,j)

            Sp1 = Energie_separare(1,1,i+1,j+1)
            Sp2 = Energie_separare(1,1,i,j)
            Sp3 = Energie_separare(1,1,i-1,j-1)

            if (Sn1[1] + Sn2[1] + Sn3[1]) == 3 # Verifica existenta simultana a celor 3 valori pentru Sₙ
                if formula == "Guttormsen"
                    Δn = abs(0.25 * (Sn1[2] - 2*Sn2[2] + Sn3[2]))
                elseif formula == "Vladuca"
                    Δn = abs(Sn2[2] - Sn3[2])
                else
                    println("==============================")
                    println("!EROARE!: Formula furnizata este gresita")
                    println(" ")
                    break
                end
                # i-p
                if isodd(j) && isodd(i)
                    push!(pairing.P, Δn)
                    push!(pairing.σᴾ, 1)
                    push!(pairing.A, i)
                    push!(pairing.Z, j)
                end
            end
            if (Sp1[1] + Sp2[1] + Sp3[1]) == 3 # Verifica existenta simultana a celor 3 valori pentru Sₚ
                if formula == "Guttormsen"
                    Δp = abs(0.25 * (Sp1[2] - 2*Sp2[2] + Sp3[2]))
                elseif formula == "Vladuca"
                    Δp = abs(Sp2[2] - Sp3[2])
                else
                    println("==============================")
                    println("!EROARE!: Formula furnizata este gresita")
                    println(" ")
                    break
                end
                # p-i
                if iseven(j) && isodd(i)
                    push!(pairing.P, Δp)
                    push!(pairing.σᴾ, 1)
                    push!(pairing.A, i)
                    push!(pairing.Z, j)
                end
            end
            # p-p
            if iseven(j) && iseven(i) && ((Sn1[1] + Sn2[1] + Sn3[1]) == 3) && ((Sp1[1] + Sp2[1] + Sp3[1]) == 3)
                push!(pairing.P, (Δn + Δp))
                push!(pairing.σᴾ, 1)
                push!(pairing.A, i)
                push!(pairing.Z, j)
            end
        end
    end
    return pairing
end
# Aici se opreste partea de calcul a programului

# Constructia reprezentarilor grafice
function Grafic_pairing(pairing, scalare, formula)
    upper_bound = maximum(pairing.P)*scalare
    if minimum(pairing.P) < 0 
        lower_bound = minimum(pairing.P)*scalare
    else
        lower_bound = 0
    end
    plt = scatter(
        pairing.A[(iseven.(pairing.A) .== 1) .& (iseven.(pairing.Z) .== 1)], 
        pairing.P[(iseven.(pairing.A) .== 1) .& (iseven.(pairing.Z) .== 1)],  
        marker = :xcross,
        markerstrokewidth = 0,
        markersize = 4, 
        ylims = (lower_bound, upper_bound),
        xlabel = L"\mathrm{A}", 
        ylabel = "Energia de împerechere [MeV]", 
        framestyle = :box,
        label = "par-par",
        title = "Energia de împerechere în funcție de A, formula $(formula)",
        minorgrid = :true,
        mc = :red, 
        size = (1000, 1000)
        )
        plt = scatter!(
        pairing.A[(iseven.(pairing.A) .== 0) .& (iseven.(pairing.Z) .== 1)], 
        pairing.P[(iseven.(pairing.A) .== 0) .& (iseven.(pairing.Z) .== 1)],  
        marker = :utriangle,
        markersize = 4, 
        markerstrokewidth = 0,
        xlabel = L"\mathrm{A}", 
        ylabel = "Energia de împerechere [MeV]", 
        framestyle = :box,
        label = "par-impar",
        minorgrid = :true,
        mc = :green, 
        size = (1000, 1000)
        )
        plt = scatter!(
        pairing.A[(iseven.(pairing.A) .== 0) .& (iseven.(pairing.Z) .== 0)], 
        pairing.P[(iseven.(pairing.A) .== 0) .& (iseven.(pairing.Z) .== 0)],  
        marker = :star5,
        markersize = 4, 
        markerstrokewidth = 0,
        xlabel = L"\mathrm{A}", 
        ylabel = "Energia de împerechere [MeV]", 
        framestyle = :box,
        label = "impar-par",
        minorgrid = :true,
        mc = :lightskyblue, 
        size = (1000, 1000)
        )     
        if lower_bound < 0
            hline!([0], ls = :dashdot, label = "")
        end
        y = 2 * 12 ./sqrt.(pairing.A)
    plt = plot!(
            pairing.A,
            y,
            xlabel = L"\mathrm{A}", 
            ylabel = "Energia de împerechere [MeV]", 
            framestyle = :box,
            label = L"2 * \frac{12}{\sqrt{A}}",
            minorgrid = :true, 
            size = (1000, 1000)
    )
    y = 12 ./sqrt.(pairing.A)
    plt = plot!(
            pairing.A,
            y,
            xlabel = L"\mathrm{A}", 
            ylabel = "Energia de împerechere [MeV]", 
            framestyle = :box,
            label = L"\frac{12}{\sqrt{A}}",
            minorgrid = :true, 
            size = (1000, 1000)
    )
    display(plt)
    #savefig(plt, "Grafice\\Pairing_$(formula).png")
end

# Apelarea functiilor definite pentru executia programului
formula = "Guttormsen"
Grafic_pairing(Energie_imperechere(formula), 1.05, formula)

formula = "Vladuca"
Grafic_pairing(Energie_imperechere(formula), 1.05, formula)