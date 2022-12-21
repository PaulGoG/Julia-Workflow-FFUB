using Plots
using CSV
using DataFrames
using LaTeXStrings

# Cod de calcul pentru partitionarea energiei totale de excitatie TXE spre fragmentele L & H

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Citire fisiere de date
df = CSV.File("Data_files/AUDI95.csv"; delim=' ', ignorerepeated=true, header=["Z", "A", "Sym", "D", "σD"]) |> DataFrame
dy = CSV.File("Data_files/U5YAZTKE.csv"; delim=' ', ignorerepeated=true, header=["A_H", "Z_H", "TKE", "Y", "σY"]) |> DataFrame

struct distributie_unidym
    x
    y
    σ
end

struct distributie_bidym
    x_1
    x_2
    y
    σ
end