#Designation of input values and parameters used throughout the main program
#####
#Add relative folder path
cd(@__DIR__)
cd("input_data/")

#=
Mass excess data file name, column names, delimiter symbol, number of first row of actual data in file
!!Header names can be rearranged, added or eliminated in order to match the actual data format of the file, but 
it should not deviate from the standard convention of eligible fields: Z, A, Symbol, D, σ_D!!
=#
mass_excess_filename = "AUDI2021.ANA"
mass_excess_header = ["Z", "A", "Symbol", "D", "σ_D"]
mass_excess_delimiter = ' '
mass_excess_firstdataline = 1

#(A,Z) of the fissionant nucleus
A₀ = 235
Z₀ = 92

#Heavy fission fragment range 
A_H_min = 118
A_H_max = 160

#Total Kinetic Energy range and step in MeV
TKE_min = 100.0
TKE_max = 200.0
TKE_step = 0.5

#=
Fission type: 
SF for spontaneous fission
(n,f) for neutron-induced fission
=#
fission_type = "(n,f)"

#Specify incident neutron energy in MeV
Eₙ = 25*1e-9

#Input variables corrections according to fission type
if fission_type == "SF"
    #Null neutron incident energy in spontaneous fission
    Eₙ = 0.0
elseif fission_type == "(n,f)"
    #Taking into account compound nucleus formation
    A₀ += 1
end

#=
Density level parameter computation method:
GC for Gilbert-Cameron
BSFG for Egidy-Bucurescu
=#
density_parameter_type = "GC"
density_parameter_filename = "SZSN.GC"

#=
Neutron evaporation cross section type:
CONSTANT for energy-independent cross section
VARIABLE for energy-dependent cross section

!If variable σ is used, provide the necessary data file for the force function S₀!
=#
evaporation_cs_type = "CONSTANT"

ħc = 197.3268601
amu = 931.50176
aₘ = 1.008665
r₀ = 1.2
global C_α = (π*ħc)^2 /(aₘ*amu)
S₀_datafile = ""

#Yield distribution filename
yield_distribution_filename = "U5YATKE.SRE"

#=
Input type for the isobaric charge distribution p(Z,A):
MEAN_VALUES for ΔZ(A)≊|0.5| & rms(A)≊0.6
DATA for values provided in a datafile

!If variable ΔZ(A) & rms(A) are used, provide the necessary data file!
=#
isobaric_distribution_type = "MEAN_VALUES"
isobaric_distribution_datafile = ""

#=
Number of Z per A fragments considered
=#
No_ZperA = 3

#=
Total Excitation Energy partitioning method:
MSCZ for Modelling at scission
PARAM for file-provided partitioning ratios

!Necessary data files must be provided in each case!: 
Extra deformation energies for MSCZ
Partitioning ratios for PARAM
=#
txe_partitioning_type = "MSCZ"
txe_partitioning_datafile = "EXTRADEF.DSE"

#Main output file name
output_filename = "dse_main_$(density_parameter_type)_$(evaporation_cs_type)_$(txe_partitioning_type).out"