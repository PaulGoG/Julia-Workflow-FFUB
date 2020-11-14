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

Sw.UAV_servo_sent_position_byte = 88 .- Sw.UAV_servo_sent_position_byte #Corectie unghi SWING

#Transformare timp fractionar in UTC

h1 = convert.(Int, last.(modf.(Md.Fractional_time)))
h2 = convert.(Int, last.(modf.(Sw.Fractional_time)))
m1 = convert.(Int, last.(modf.(first.(modf.(Md.Fractional_time)) .*60)))
m2 = convert.(Int, last.(modf.(first.(modf.(Sw.Fractional_time)) .*60)))
s1 = convert.(Int, round.((first.(modf.( first.(modf.(Md.Fractional_time)) .*60) )) .* 60))
s2 = convert.(Int, round.((first.(modf.( first.(modf.(Sw.Fractional_time)) .*60) )) .* 60))

timp = Timp(string.(h1, ':', m1, ':', s1), string.(h2, ':', m2, ':', s2))

unghiuri_Md = unique(Md.Elev_viewing_angle)
unghiuri_Sw = sort(unique(Sw.UAV_servo_sent_position_byte))

#Plot MaxDoas toate unghiurile

plt = scatter(framestyle =:box,
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
    plt = scatter!(x, y, label = "α = $(unghiuri_Md[i])ᵒ",
    size=(1280,900))
end

display(plt)
savefig("MaxDoas.pdf")

#Plot SWING toate unghiurile

plt = scatter(framestyle =:box,
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
    plt = scatter!(x, y, label = "α = $(unghiuri_Sw[i])ᵒ",
    size=(1280,900))
end

display(plt)
savefig("Swing.pdf")

#Plot ambele dispozitive la unghiurile comune

for i in 1:length(unghiuri_Md)
    j = findfirst(x -> x == unghiuri_Md[i], unghiuri_Sw)

    if j !== nothing
        
    indici_Md = findall(x -> x == unghiuri_Md[i], Md.Elev_viewing_angle)
    indici_Sw = findall(x -> x == unghiuri_Sw[j], Sw.UAV_servo_sent_position_byte)
    
    x = Md.Fractional_time[indici_Md]
    y = Md.NO2_SlCol_no2_[indici_Md]
    plt = scatter(x, y, label = "MaxDoas")

    x = Sw.Fractional_time[indici_Sw]
    y = Sw.NO2_SlCol_NO2_[indici_Sw]
    plt = scatter!(x, y, label = "Swing")

    plt = scatter!(framestyle =:box,
    title = "%c NO₂ - la unghiul α = $(unghiuri_Md[i])ᵒ",
    legend = :outertopright,
    xlabel = "Timp (UTC)",
    ylabel = "Concentratie NO₂",
    xticks = ( Md.Fractional_time[collect(1:49:length(Md.Fractional_time))], timp.Md[collect(1:49:length(Md.Fractional_time))] ),
    xrotation = 60,
    size=(1280,900)
    );

    display(plt)
    savefig("Ambele_dispozitive_$(unghiuri_Md[i]).pdf")
    end
end
