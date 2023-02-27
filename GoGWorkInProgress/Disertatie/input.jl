#Designation of input values and parameters used throughout the main program
#####
#Add path to input data folder
cd(@__DIR__)
cd("input_data/")

#Mass excess data file name, column names, delimiter symbol, number of first row of actual data in file

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
TKE_step = 2.0

#=
Fission type: 
SF for spontaneous fission
(n,f) for neutron-induced fission
=#
fission_type = "(n,f)"

#Specify incident neutron energy in MeV
Eₙ = 25*1e-9

#=
Density level parameter computation method:
GC for Gilbert-Cameron
BSFG for Egidy-Bucurescu
=#
density_parameter_type = "BSFG"
density_parameter_filename = "SZSN.GC"
density_parameter_header = ["n", "S_Z", "S_N"]
density_parameter_delimiter = ' '
density_parameter_firstdataline = 2

#=
Neutron evaporation cross section type:
CONSTANT for energy-independent cross section
VARIABLE for energy-dependent cross section
=#
evaporation_cs_type = "VARIABLE"

#=
Input type for the isobaric charge distribution p(Z,A):
MEAN_VALUES for ΔZ(A)=0.5 & rms(A)=0.6
DATA for values provided in a datafile

!If variable ΔZ(A) & rms(A) are used, provide the necessary data file!
=#
isobaric_distribution_type = "DATA"
isobaric_distribution_filename = "DeltaZA_rmsA.U5"
isobaric_distribution_header = ["A", "ΔZ_A", "rms_A"]
isobaric_distribution_delimiter = ' '
isobaric_distribution_firstdataline = 2

#=
Number of Z per A fragments considered
=#
No_ZperA = 5

#=
Total Excitation Energy partitioning method:
MSCZ for Modelling at scission
PARAM for file-provided partitioning ratios

!Necessary data files must be provided in each case!: 
Extra deformation energies for MSCZ
Partitioning ratios for PARAM
=#
txe_partitioning_type = "MSCZ"
txe_partitioning_filename = "EXTRADEF.IN"
txe_partitioning_header = ["A", "Z", "Value"]
txe_partitioning_delimiter = ' '
txe_partitioning_firstdataline = 2

#Main output file name
output_filename = "DSE_$(No_ZperA)ZA_$(isobaric_distribution_type)_pAZ_$(density_parameter_type)_cs_$(evaporation_cs_type)_partitioning_$(txe_partitioning_type).OUT"