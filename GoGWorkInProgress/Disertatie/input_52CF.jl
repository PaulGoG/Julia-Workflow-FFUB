#Designation of input values and parameters 
#####
using Dates
println("*begin program execution at $(Dates.format(now(), "HH:MM:SS"))")
println("*loading input parameters and Julia libraries")
#Add path to input data folder
cd(@__DIR__)
if !isdir("input_data/")
    mkdir("input_data/")
end
cd("input_data/")

#Mass excess data file name, column names, delimiter symbol, number of first row of actual data in file
mass_excess_filename = "AME2020.ANA"
mass_excess_header = ["Z", "A", "Symbol", "D", "σ_D"]
mass_excess_delimiter = ' '
mass_excess_firstdataline = 1

#(A,Z) and symbol identifier of the fissionant nucleus
fissionant_nucleus_identifier = "CF52"
A₀ = 252
Z₀ = 98

#Heavy fission fragment range 
A_H_min = 126
A_H_max = 174

#Total Kinetic Energy range and step in MeV
TKE_min = 130.0
TKE_max = 230.0
TKE_step = 2.0

#=
Fission type: 
*SF for spontaneous fission
*(n,f) for neutron-induced fission
=#
fission_type = "SF"

#Specify incident particle energy in MeV
E_incident = 0.0

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
*VARIABLE for energy-dependent cross section modelled with s-wave neutron strength function
=#
evaporation_cs_type = "VARIABLE"

#=
Input type for the isobaric charge distribution p(Z,A):
*MEAN_VALUES for ΔZ(A_H) = -0.5 & rms(A_H) = 0.6
*DATA for values provided in a datafile
=#
isobaric_distribution_type = "DATA"

isobaric_distribution_filename = "DeltaZ_rms_A.$fissionant_nucleus_identifier"
isobaric_distribution_header = ["A", "ΔZ_A", "rms_A"]
isobaric_distribution_delimiter = ' '
isobaric_distribution_firstdataline = 2

#Number of normally distributed Z considered for each A
No_ZperA = 5

#=
Total Excitation Energy partitioning method:
*MSCZ for Modelling at scission
*PARAM for file-provided partitioning ratios E*_H/TXE via segments as a function of A_H 
*RT as a function of A_H for T_L/T_H ratio provided by user via segments -for constant RT provide y=const. segment line-

Data provided in each case: 
*Extra deformation energies for MSCZ via datafile
*Segment points in Vector of Tuples (txe_partitioning_segmentpoints) for PARAM & RT
=#
txe_partitioning_type = "PARAM"

txe_partitioning_filename = ""
txe_partitioning_header = ["A", "Z", "Value"]
txe_partitioning_delimiter = ' '
txe_partitioning_firstdataline = 2

txe_partitioning_segmentpoints = [(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]

#Yield-averaged outputs YES or NO selector for data and plots
secondary_outputs = "YES"

secondary_output_Yield = "YES"
secondary_output_TXE_Q = "YES"
secondary_output_E_excitation = "YES"
secondary_output_ν = "YES"
secondary_output_Ap = "YES"
secondary_output_T = "YES"
secondary_output_Tₖ = "YES"
secondary_output_avg_ε = "NO"
secondary_output_avg_εₖ = "NO"
secondary_output_Eᵣ = "NO"
ΔT, ΔTₖ, Δavg_ε, Δavg_εₖ, ΔEᵣ = 1e-2, 1e-2, 5e-2, 5e-2, 5e-1    #Histogram steps

yield_distribution_filename = "$(fissionant_nucleus_identifier)YATKE.VES"
yield_distribution_header = ["A", "TKE", "Value", "σ"]
yield_distribution_delimiter = ' '
yield_distribution_firstdataline = 2

#Neutron spectrum calculation YES or NO selector   
neutron_spectrum = "NO"
E_min = 1e-6
E_max = 20.0
E_step = 5e-2
#Yield_cutoff_value = 

#Plots YES or NO selector
generate_plots = "NO"
multidim_plots = "NO"

resolution_scale = 250
aspect_ratio = (3, 2)