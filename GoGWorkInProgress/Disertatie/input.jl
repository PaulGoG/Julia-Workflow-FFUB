#Designation of input values and parameters used throughout the main program
#####
#Add relative folder path
cd(@__DIR__);

#Mass excess data file name and relative path
mass_excess_filename = "input_data/AUDI2021.ANA";

#(A,Z) of the fissionant nucleus
A₀,Z₀ = 92,236;

#Heavy fission fragment range 
A_H_min,A_H_max = 118,160;

#Total Kinetic Energy range and step in MeV
TKE_min,TKE_max,TKE_step = 100.0,200.0,0.5;

#=
Fission type: 
SF for spontaneous fission
(n,f) for neutron-induced fission
=#
fission_type = "(n,f)";
if fission_type == "SF";
    #Null neutron incident energy in spontaneous fission
    Eₙ = 0.0;
elseif fission_type == "(n,f)";
    #Specify incident neutron energy in MeV
    Eₙ = 25*1e-9;
    #Taking into account compound nucleus formation
    A₀ += 1;
else println("=====!ERROR, WRONG INPUT FOR fission_type VARIABLE!=====");
end

#=
Density level parameter computation method:
GC for Gilbert-Cameron
BSFG for Egidy-Bucurescu

!Necessary data files must be provided in each case!
=#
density_parameter_type = "GC";
if density_parameter_type == "GC";
    density_parameter_filename = "input_data/SZSN.GC";
elseif density_parameter_type == "BSFG";
    density_parameter_filename = "input_data/";
else println("=====!ERROR, WRONG INPUT FOR density_parameter_type VARIABLE!=====");
end

#=
Neutron evaporation cross section type:
CONSTANT for energy-independent cross section
VARIABLE for energy-dependent cross section

!If variable σ is used, provide the necessary data file for the force function S₀!
=#
evaporation_cs_type = "CONSTANT";
if evaporation_cs_type == "VARIABLE";
    const ħc = 197.3268601;
    const a.m.u = 931.50176;
    const aₘ = 1.008665;
    const r₀ = 1.2;
    const C_α = (π*ħc)^2 /(aₘ*a.m.u);
    S₀_datafile = "input_data/";
end

#Yield distribution path and filename
yield_distribution_filename = "input_data/U5YATKE.SRE";

#=
Input type for the isobaric charge distribution p(Z,A):
MEAN_VALUES for ΔZ(A)≊|0.5| & rms(A)≊0.6
DATA for values provided in a datafile

!If variable ΔZ(A) & rms(A) are used, provide the necessary data file!
=#
isobaric_distribution_type = "MEAN_VALUES";
if isobaric_distribution_type == "DATA";
    isobaric_distribution_datafile = "input_data/";
end

#=
Number of Z per A fragments considered
!!MUST BE AN ODD INTEGER VALUE!!
=#
No_ZperA = 7;

#Main output file name
output_filename = "dse_main.out";
