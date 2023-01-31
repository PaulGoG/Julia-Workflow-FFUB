using Plots
using CSV
using DataFrames
using LaTeXStrings
using QuadGK
using Optim
using Trapz

# Cod de calcul pentru fitul datelor experimentale ale spectrelor neutronice la un spectru Maxwell

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru

# Citire fisiere de date
dn_Gook_SL = CSV.read("Data_files/Date_experimentale/Spectru_n/U5SPGOOK.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["E", "N", "σN"], skipto = 2);
dn_Vorobyev_SL = CSV.read("Data_files/Date_experimentale/Spectru_n/U5SPVORO.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["E", "N", "σN"], skipto = 2);
dn_Gook_SCM_LF = CSV.read("Data_files/Date_experimentale/Spectru_n/U5SPCMLF.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["E", "N", "σN"], skipto = 2);
dn_Gook_SCM_HF = CSV.read("Data_files/Date_experimentale/Spectru_n/U5SPCMHF.DAT", DataFrame; delim=' ', ignorerepeated=true, header=["E", "N", "σN"], skipto = 2);
#####
# Nₘ(E)
function N_Maxwell(E, T_M)
    return (2/sqrt(π)) * T_M^(-3/2) * sqrt(E) * exp(-E/T_M)
end
function Arie_date_experimentale(dn)
    return trapz(dn.E, dn.N)
end
function Arie_teoretica(E_min, E_max, T_M)
    Maxwell(E) = (2/sqrt(π)) * T_M^(-3/2) * sqrt(E) * exp(-E/T_M)
    return quadgk(Maxwell, E_min, E_max)[1]
end
function Normare_Maxwell(dn, T_M)
    E_min = first(dn.E)
    E_max = last(dn.E)
    return Arie_teoretica(E_min, E_max, T_M)/Arie_date_experimentale(dn)
end
function Normare_Spectru(dn, T_M)
    f = Normare_Maxwell(dn, T_M)
    for i in eachindex(dn.E)
        dn.N[i] *= f/N_Maxwell(dn.E[i], T_M)
        dn.σN[i] *= f/N_Maxwell(dn.E[i], T_M)
    end
    return dn
end

dn = dn_Gook_SL;
χ²(T_M) = (1/length(dn.E)) * sum([((dn.N[i] - (N_Maxwell(dn.E[i], T_M)/Normare_Maxwell(dn, T_M)))/dn.σN[i])^2 for i in eachindex(dn.N)]);
result = optimize(χ², 0.0, 3.0);
T_M = Optim.minimizer(result);
χ_sq = Optim.minimum(result);
scatter(dn.E, dn.N, yerror = dn.σN, yscale = :log10);
plot!(dn.E, N_Maxwell.(dn.E, T_M))
dn = Normare_Spectru(dn, T_M);
scatter(dn.E, dn.N, yerror = dn.σN, xscale = :log10);
hline!([1.0])