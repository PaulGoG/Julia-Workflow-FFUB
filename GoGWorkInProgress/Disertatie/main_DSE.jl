#=
Main body of the DSE (Deterministic Sequential Emission) program;
It only outputs primary results obtained by solving the sequential emission energy conservation equations;
For secondary results comparable to experimental data (eg. neutron multiplicity) run the
secondary_output.jl program, and for graphical representations run the graphics_plots.jl program;
=#
#####
#Include auxiliary files in main program
include("input.jl")
include("error_handling.jl")
include("aux_func.jl")
include("density_parameters.jl")
include("txe_partitioning.jl")

#Read input data files as DataFrames
dmass_excess = CSV.read(mass_excess_filename, DataFrame; delim = mass_excess_delimiter, ignorerepeated = true, header = mass_excess_header, skipto = mass_excess_firstdataline)

#Revert relative PATH to project root folder
cd(@__DIR__)
