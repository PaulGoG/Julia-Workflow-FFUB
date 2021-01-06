using Plots
using CSV
using DataFrames

cd(@__DIR__) #Adauga path-ul la fisier
 # Z|A|Simbol|D(KeV)|σᴰ(KeV) -> forma tabelului
df = CSV.read("AUDI95.csv"; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σ"]) #citim fisierul de tip CSV

Dn = df.D[(df.Z .== 0) .& (df.A .== 1)][1] #Defectul de masa al neutronului
D236 = df.D[(df.Z .== 92) .& (df.A .== 236)][1] #Defectul de masa al U236 nucleu compus
D235 = df.D[(df.Z .== 92) .& (df.A .== 235)][1] #Defectul de masa al U235
σn = df.σ[(df.Z .== 0) .& (df.A .== 1)][1]      #Abaterile/incertitudinile pentru fiecare
σ236 = df.σ[(df.Z .== 92) .& (df.A .== 236)][1]
σ235 = df.σ[(df.Z .== 92) .& (df.A .== 235)][1]


Qc = zeros(7) #Creem vectori goi pe care ii populam cu valori de calcul
Qnc = zeros(7)
σQc = zeros(7)
σQnc = zeros(7)

ok = 53
clearconsole()
for i in 137:143
    print("Fragmentului greu: ")
    print(df.Sym[(df.Z .== ok) .& (df.A .== i), :][1])
    print(" cu A = ", i," & Z = ", ok, "\n")
    print("Ii corespunde fragmentul usor: ")
    print(df.Sym[(df.Z .== 92 - ok) .& (df.A .== 236 - i), :][1])
    print(" cu A = ", 236-i," & Z = ", 92-ok, "\n")

    DH = df.D[(df.Z .== ok) .& (df.A .== i), :][1] #Defectul de masa al fragmentului greu
    DL = df.D[(df.Z .== 92 - ok) .& (df.A .== 236 - i), :][1] #Defectul de masa al fragmentului usor
    σH = df.σ[(df.Z .== ok) .& (df.A .== i), :][1]            #Abaterile corespunzatoare
    σL = df.σ[(df.Z .== 92 - ok) .& (df.A .== 236 - i), :][1]

    Q = (D235 + Dn - (DH + DL))/1000
    Qnc[i-136] = Q
    σQ = (σ235 + σn + σH + σL)/1000
    σQnc[i-136] = σQ
    print("Energia eliberata la fisiune fara a tine cont de formarea nucleului compus Q = " ,Q ," MeV, si abaterea/incertitudinea corespunzatoare σ = ",σQ ," MeV", "\n")
                                #Tinem cont de corectia formarii nucleului compus
    Q = (D236 - (DH + DL))/1000
    Qc[i-136] = Q
    σQ = (σ236 + σH + σL)/1000
    σQc[i-136] = σQ
    print("Energia eliberata la fisiune tinand cont de formarea nucleului compus Q = " ,Q ," MeV, si abaterea/incertitudinea corespunzatoare σ = ",σQ ," MeV", "\n", "\n")

    if (i == 138 || i == 141)
        global ok+=1
    end
end

Sn = (D235-D236+Dn)/1000 #Corectia luarii in calcul a formarii nucleului compus
println()
println("Energia pentru a scoate neutronul din nucleul compus este: ", Sn, " MeVi")


Y = [3.3504, 3.0624, 4.0473, 4.3027, 3.5286, 3.5328, 3.0646]
Y = Y./1000
σY = [0.0288, 0.0257, 0.0336, 0.0353, 0.0291, 0.0297, 0.0260]
σY = σY./1000       #Distributia dupa masa si abaterea ei in MeV

Qcmed = sum(Qc.*Y)/sum(Y)
σcmed = sum(Y.*σQc .+ Qc.*σY)/sum(Y) + sum(Qc.*Y)*sum(σY)/sum(Y)^2
println("In cazul luarii in considerare a formarii nucleului compus, <Q> = ", Qcmed, " +- ", σcmed, " MeV")
Qncmed = sum(Qnc.*Y)/sum(Y)
σncmed = sum(Y.*σQnc .+ Qnc.*σY)/sum(Y) + sum(Qnc.*Y)*sum(σY)/sum(Y)^2
println("In cazul neluarii in considerare a formarii nucleului compus, <Q> = ", Qncmed, " +- ", σncmed, " MeV")

#Corectat inseamna cu nucleu compus
#Necorectat inseamna fara nucleu compus
x = 137:143
plot(x, Qc, yerror = σQc, framestyle =:box, title = "Diferenta luarii in calcul a formarii nucleului compus", label = "QNucleuCompus", foreground_color_legend = nothing, background_color_legend = nothing)
plot!(x,Qnc, yerror = σQnc, label = "QFaraNucleuCompus")
xlabel!("A")
ylabel!("Energia eliberata la fisiune (MeV)")
y = fill(Qcmed, 7)
err = fill(σcmed, 7)
plot!(x, y, label = "QNCmediu", yerror = err)
y = fill(Qncmed, 7)
err = fill(σncmed, 7)
plot!(x, y, label = "QFNCmediu", yerror = err)


# sqrt(sum( ((Y./sum(Y)).^2).*(σQc.^2)) + sum((((Qc.-Qcmed)./sum(Y)).^2).*(σY.^2)))
# Formula Propagarii erorilor calculata de profa cu derivate + ridicare la patrat
# σ trebuie corectate cu formula abaterii patratice medii !!! (For future use)
