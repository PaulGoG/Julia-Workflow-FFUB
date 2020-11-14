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

#Ploturi

unghiuri_Md = unique(Md.Elev_viewing_angle)
unghiuri_Sw = sort(unique(Sw.UAV_servo_sent_position_byte))

scatter(framestyle =:box,
title = "MaxDoas",
foreground_color_legend = nothing, background_color_legend = nothing)

for i in 1:length(unghiuri_Md)
    indici = findall(x -> x == unghiuri_Md[i], Md.Elev_viewing_angle)
    x = Md.Fractional_time[indici]
    y = Md.NO2_SlCol_no2_[indici]
    scatter(x, y, label = "α = $(unghiuri_Md[i])")
end

    indici = findall(x -> x == unghiuri_Md[6], Md.Elev_viewing_angle)
    x = Md.Fractional_time[indici]
    y = Md.NO2_SlCol_no2_[indici]
    scatter(x, y, label = "α = $(unghiuri_Md[6])")

    unghiuri_Md[8]


scatter(Md.Fractional_time, Md.NO2_SlCol_no2_ , framestyle =:box,
title = "MaxDoas", label = "Label",
foreground_color_legend = nothing, background_color_legend = nothing)

plot!(x,Qnc, yerror = σQnc, label = "QFaraNucleuCompus")
xlabel!("A")
ylabel!("Energia eliberata la fisiune (MeV)")
y = fill(Qcmed, 7)
err = fill(σcmed, 7)
plot!(x, y, label = "QNCmediu", yerror = err)
y = fill(Qncmed, 7)
err = fill(σncmed, 7)
plot!(x, y, label = "QFNCmediu", yerror = err)
