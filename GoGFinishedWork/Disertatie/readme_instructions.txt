Instructions for input.jl and options.jl file parameters
########################################################################################################

!This public version of the program has been redacted and is missing key files and function definitions necessary for it to run succesfully!

########################################################################################################

Variables used for automating data reading from files have the following structure:
data file name, column names, delimiter symbol, number of first row of actual data in file

Fission type options: 
*SF for spontaneous fission
*(n,f) for neutron-induced fission

Density level parameter computation method:
*GC for Gilbert-Cameron using shell corrections
*BSFG for Egidy-Bucurescu Back Shifted Fermi Gas

Neutron evaporation cross section type:
*CONSTANT for constant cross section
*VARIABLE for energy-dependent cross section modelled with s-wave neutron strength function

Input type for the isobaric charge distribution p(Z,A):
*MEAN_VALUES for ΔZ(A_H) = -0.5 & rms(A_H) = 0.6
*DATA for values provided in a data file

Total Excitation Energy (TXE) partitioning method:
*MSCZ for Modelling at scission
*PARAM for file-provided partitioning ratios E*_H/TXE via segments as a function of A_H 
*RT as a function of A_H for T_L/T_H ratio provided by user via segments -for constant RT provide y=const. segment line-
Data provided in each case: 
*Extra deformation energies for MSCZ via data file
*Segment points in Vector of Tuples (txe_partitioning_segmentpoints) for PARAM & RT