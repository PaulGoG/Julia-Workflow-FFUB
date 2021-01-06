
#Pachetele trebuie instalate din consola la prima rulare cu comanda Pkg.add !!!

using Plots
using CSV
using DataFrames

mutable struct Timp
    Md
    Sw
end

cd(@__DIR__) #Adauga path-ul relativ la folder

Md = CSV.File("MAXDOAS.csv"; select=[5, 8, 11], normalizenames=true) |> DataFrame  #Citim coloanele de
Sw = CSV.File("SWING.csv"; select=[3, 6, 14], normalizenames=true) |> DataFrame    #interes din fisiere

Sw.UAV_servo_sent_position_byte = 88 .- Sw.UAV_servo_sent_position_byte #Corectie unghi SWING (operatii vectoriale)

#Transformare timp fractionar in UTC
#Pentru timp fractionar denumit T.frac
# h = [T.frac]
# m = [{T.frac} * 60]
# s = {{T.frac} * 60} * 60

h1 = string.(convert.(Int, last.(modf.(Md.Fractional_time))))
h2 = string.(convert.(Int, last.(modf.(Sw.Fractional_time))))
m1 = string.(convert.(Int, last.(modf.(first.(modf.(Md.Fractional_time)) .*60))))
m2 = string.(convert.(Int, last.(modf.(first.(modf.(Sw.Fractional_time)) .*60))))
s1 = string.(convert.(Int, round.((first.(modf.( first.(modf.(Md.Fractional_time)) .*60) )) .* 60)))
s2 = string.(convert.(Int, round.((first.(modf.( first.(modf.(Sw.Fractional_time)) .*60) )) .* 60)))

#Adaugam 0-uri in fata cifrelor

for i in 1:length(Md.Fractional_time)
    if parse(Int, h1[i]) < 10
        h1[i] = lpad(h1[i], 2, "0")
    end
    if parse(Int, m1[i]) < 10
        m1[i] = lpad(m1[i], 2, "0")
    end
    if parse(Int, s1[i]) < 10
        s1[i] = lpad(s1[i], 2, "0")
    end
end

for i in 1:length(Sw.Fractional_time)
    if parse(Int, h2[i]) < 10
        h2[i] = lpad(h2[i], 2, "0")
    end
    if parse(Int, m2[i]) < 10
        m2[i] = lpad(m2[i], 2, "0")
    end
    if parse(Int, s2[i]) < 10
        s2[i] = lpad(s2[i], 2, "0")
    end
end

#Creem obiectul de tip struct cu timpii in format UTC

timp = Timp(string.(h1, ':', m1, ':', s1), string.(h2, ':', m2, ':', s2))

unghiuri_Md = unique(Md.Elev_viewing_angle)
unghiuri_Sw = sort(unique(Sw.UAV_servo_sent_position_byte))

#Grafic MaxDoas cu toate unghiurile

plt = plot(
    framestyle =:box,
    title = "%c NO₂ - MaxDoas",
    legend = :outertopright,
    xlabel = "Timp (UTC)",
    ylabel = "Concentratie NO₂",
    xticks = ( Md.Fractional_time[collect(1:50:length(Md.Fractional_time))], timp.Md[collect(1:50:length(Md.Fractional_time))] )
);

for i in 1:length(unghiuri_Md)
    indici = findall(x -> x == unghiuri_Md[i], Md.Elev_viewing_angle)
    x = Md.Fractional_time[indici]
    y = Md.NO2_SlCol_no2_[indici]
    plt = plot!(x, y, label = "α = $(unghiuri_Md[i])ᵒ",
    size = (1280,900), shape = :auto)
end

display(plt)
savefig("Grafice\\MaxDoas.png")

#Grafic SWING toate unghiurile

plt = plot(
    framestyle =:box,
    title = "%c NO₂ - Swing",
    legend = :outertopright,
    xlabel = "Timp (UTC)",
    ylabel = "Concentratie NO₂",
    xticks = ( Sw.Fractional_time[collect(1:70:length(Sw.Fractional_time))], timp.Sw[collect(1:70:length(Sw.Fractional_time))] )
);

for i in 1:length(unghiuri_Sw)
    indici = findall(x -> x == unghiuri_Sw[i], Sw.UAV_servo_sent_position_byte)
    x = Sw.Fractional_time[indici]
    y = Sw.NO2_SlCol_NO2_[indici]
    plt = plot!(x, y, label = "α = $(unghiuri_Sw[i])ᵒ",
    size = (1280,900), shape = :auto)
end

display(plt)
savefig("Grafice\\Swing.png")

#Grafice cu ambele dispozitive la toate unghiurile comune

for i in 1:length(unghiuri_Md)
    j = findfirst(x -> x == unghiuri_Md[i], unghiuri_Sw) #Gasim unghiurile comune

    if j !== nothing #IMPORTANT!!!
        
    indici_Md = findall(x -> x == unghiuri_Md[i], Md.Elev_viewing_angle)
    indici_Sw = findall(x -> x == unghiuri_Sw[j], Sw.UAV_servo_sent_position_byte)
    
    x = Md.Fractional_time[indici_Md]
    y = Md.NO2_SlCol_no2_[indici_Md]
    plt = plot(x, y, label = "MaxDoas", shape = :auto)

    x = Sw.Fractional_time[indici_Sw]
    y = Sw.NO2_SlCol_NO2_[indici_Sw]
    plt = plot!(x, y, label = "Swing", shape = :auto)

    plt = plot!(
        framestyle =:box,
        title = "%c NO₂  la unghiul α = $(unghiuri_Md[i])ᵒ",
        legend = :outertopright,
        xlabel = "Timp (UTC)",
        ylabel = "Concentratie NO₂",
        xticks = ( Md.Fractional_time[collect(1:49:length(Md.Fractional_time))], timp.Md[collect(1:49:length(Md.Fractional_time))] ),
        xrotation = 60,
        size = (1280, 900),
    );

    display(plt)
    savefig("Grafice\\Ambele_dispozitive_$(unghiuri_Md[i]).png")
    end
end