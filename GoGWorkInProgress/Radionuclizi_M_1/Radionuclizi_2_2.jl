using CSV
using DataFrames

cd(@__DIR__) #Adauga path-ul la fisier

function Calcul_citire_tastatura()
    df = CSV.File("AUDI95.csv"; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) |> DataFrame
    global contor = true
    while(contor)
        global contor_psep = true
        global contor_nuc = true
        while(contor_psep)
            # Citire A_particula si Z_particula
            println(" ")
            println("======================================================================================================================")
            println("Introduceți valorile corespunzătoare pentru A-ul și Z-ul particulei ce urmează a fi separată din nucleu: ")
            println(" ")
            print("A_particulă = ")
            global A_part = parse(Int, readline())
            print("Z_particulă = ")
            global Z_part = parse(Int, readline())
            if isassigned(df.D[(df.A .== A_part) .& (df.Z .== Z_part)], 1)
                global contor_psep = false
                println(" ")
                break
            else 
                println(" ")
                println(" ")
                println(" ")
                println("======================================================!ATENȚIE!=======================================================")
                println("Particula cu (A = $(A_part), Z = $(Z_part)) introduse de la tastatură nu există în librăria de date nucleare!")
                println("Urmează reluarea procesului de introducere a datelor de la tastatură!")
                println("======================================================================================================================")
                println(" ")
                println(" ")
                println(" ")
            end
        end
        while(contor_nuc)
            # Citire A si Z ale nucleului din care vrem sa separam particula
            println("Introduceți valorile corespunzătoare pentru A-ul și Z-ul nucleului din care separăm particula $(df.Sym[(df.A .== A_part) .& (df.Z .== Z_part)][1])($(A_part), $(Z_part)):")
            println(" ")
            print("A = ")
            global A = parse(Int, readline())
            print("Z = ")
            global Z = parse(Int, readline())
            if isassigned(df.D[(df.A .== A) .& (df.Z .== Z)], 1)
                global contor_nuc = false
                println(" ")
            else 
                println(" ")
                println(" ")
                println(" ")
                println("======================================================!ATENȚIE!=======================================================")
                println("Nucleul cu (A = $A, Z = $Z) introduse de la tastatură nu există în librăria de date nucleare!")
                println("Urmează reluarea procesului de introducere a datelor de la tastatură pentru nucleul din care separăm particula $(df.Sym[(df.A .== A_part) .& (df.Z .== Z_part)][1])($(A_part), $(Z_part))!")
                println("======================================================================================================================")
                println(" ")
                println(" ")
                println(" ")
            end

            if contor_nuc == false
                if isassigned(df.D[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)], 1)
                    global contor = false
                    break
                else
                    println(" ")
                    println(" ")
                    println(" ")
                    println("======================================================!ATENȚIE!=======================================================")
                    println("Nucleul rezultat în urma separării având (A = $(A - A_part), Z = $(Z - Z_part)) nu există în librăria de date nucleare!")
                    println("Urmează reluarea procesului de introducere a datelor de la tastatură!")
                    println("======================================================================================================================")
                    println(" ")
                    println(" ")
                    println(" ")
                end
            end
        end
    end

    D_part = df.D[(df.A .== A_part) .& (df.Z .== Z_part)][1]
    σ_part = df.σ[(df.A .== A_part) .& (df.Z .== Z_part)][1]

    D = df.D[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)][1]
    σᴰ = df.σ[(df.A .== A - A_part) .& (df.Z .== Z - Z_part)][1]    

    S = round((D + D_part - df.D[(df.A .== A) .& (df.Z .== Z)][1])*1e-3, digits = 7)
    σˢ = round(sqrt(σ_part^2 + σᴰ^2 + (df.σ[(df.A .== A) .& (df.Z .== Z)][1])^2)*1e-3, digits = 7)

    println("Energia de separare a $(df.Sym[(df.A .== A_part) .& (df.Z .== Z_part)][1])($(A_part), $(Z_part)) din nucleul de $(df.Sym[(df.A .== A) .& (df.Z .== Z)][1])($(A), $(Z)) este de $(S) +- $(σˢ) [MeV]")
    println("======================================================================================================================")
    println(" ")
end

Calcul_citire_tastatura()