# Aici citim CSV-urile cu date tabelate cu care o sa efectuam calculele

cd(@__DIR__)

T_1 = CSV.File("Date_Tabelate_CSV\\Tabel_1.csv"; normalizenames=true) |> DataFrame
T_2 = CSV.File("Date_Tabelate_CSV\\Tabel_2.csv"; normalizenames=true) |> DataFrame