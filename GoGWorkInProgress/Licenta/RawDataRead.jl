# Aici citim CSV-urile cu date tabelate cu care o sa efectuam calculele

cd(@__DIR__)

T_1 = CSV.File("Date_Tabelate_CSV\\Tabel_1.csv"; normalizenames=true) |> DataFrame
T_2 = CSV.File("Date_Tabelate_CSV\\Tabel_2.csv"; normalizenames=true) |> DataFrame
T_3 = CSV.File("Date_Tabelate_CSV\\Tabel_3.csv"; normalizenames=true) |> DataFrame
T_4 = CSV.File("Date_Tabelate_CSV\\Tabel_4.csv"; normalizenames=true) |> DataFrame
T_7 = CSV.File("Date_Tabelate_CSV\\Tabel_7.csv"; normalizenames=true) |> DataFrame
Cladiri = CSV.File("Date_Tabelate_CSV\\Cladiri.csv"; normalizenames=true) |> DataFrame
Freq = CSV.File("Date_Tabelate_CSV\\Frecvente.csv"; normalizenames=true) |> DataFrame