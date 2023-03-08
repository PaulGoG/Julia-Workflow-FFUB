#Designation of input values and parameters 
#####
#Add path to input data folder
cd(@__DIR__)
if !isdir("input_data/")
    mkdir("input_data/")
end
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
*GC for Gilbert-Cameron using shell corrections
*BSFG for Egidy-Bucurescu Back Shift Fermi Gas
=#
density_parameter_type = "BSFG"

density_parameter_filename = "SZSN.GC"
density_parameter_header = ["n", "S_Z", "S_N"]
density_parameter_delimiter = ' '
density_parameter_firstdataline = 2

#=
Neutron evaporation cross section type:
*CONSTANT for constant cross section
*VARIABLE for energy-dependent cross section modelled with s-wave neutron force function
=#
evaporation_cs_type = "VARIABLE"

#=
Input type for the isobaric charge distribution p(Z,A):
*MEAN_VALUES for ΔZ(A_H) = -0.5 & rms(A_H) = 0.6
*DATA for values provided in a datafile
=#
isobaric_distribution_type = "DATA"
isobaric_distribution_filename = "DeltaZA_rmsA.U5"
isobaric_distribution_header = ["A", "ΔZ_A", "rms_A"]
isobaric_distribution_delimiter = ' '
isobaric_distribution_firstdataline = 2

#Number of normally distributed Z considered for each A
No_ZperA = 5

#=
Total Excitation Energy partitioning method:
*MSCZ for Modelling at scission
*PARAM for file-provided partitioning ratios E*_H/TXE via segments
*RT(A_H) for T_L/T_H ratio provided by user via segments -for constant RT provide segment line-

!Data must be provided in each case!: 
*Extra deformation energies for MSCZ via datafile
*Segment points in Vector of Tuples for PARAM & RT
=#
txe_partitioning_type = "MSCZ"
txe_partitioning_filename = "EXTRADEF.IN"
txe_partitioning_header = ["A", "Z", "Value"]
txe_partitioning_delimiter = ' '
txe_partitioning_firstdataline = 2
txe_partitioning_segmentpoints = [(118, 1.2), (160, 1.2)]

#Writing out main DSE output file containing detailed sequence data YES or NO selector
write_primary_output = "NO"

#Yield-averaged outputs YES or NO selector
secondary_outputs = "YES"
yield_distribution_filename = "U5YATKE.ANA"
yield_distribution_header = ["A", "TKE", "Value", "σ"]
yield_distribution_delimiter = ' '
yield_distribution_firstdataline = 2

#=
Plots YES or NO selector
!It requires YES to secondary_outputs!
=#
generate_plots = "NO"

#Neutron spectrum calculation YES or NO selector   
neutron_spectrum = "NO"