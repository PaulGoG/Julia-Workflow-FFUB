#=
Main body of the DSE (Deterministic Sequential Emission) program;
It only outputs primary results obtained by solving the sequential emission energy conservation equations;
For secondary results comparable to experimental data (eg. neutron multiplicity) run the
secondary_output.jl program, and for graphical representations run the graphics_plots.jl program;
=#
#####
println("*starting program")
println("*loading input parameters and Julia libraries")
#Include auxiliary files in main program
include("input.jl")
include("error_handling.jl")
include("aux_func.jl")
include("density_parameters.jl")
include("txe_partitioning.jl")

println("*reading data files")
#Read input data files as DataFrames
dmass_excess = CSV.read(mass_excess_filename, DataFrame; delim = mass_excess_delimiter, ignorerepeated = true, header = mass_excess_header, skipto = mass_excess_firstdataline)

if density_parameter_type == "GC"
    dGC = CSV.read(density_parameter_filename, DataFrame; delim = density_parameter_delimiter, ignorerepeated = true, header = density_parameter_header, skipto = density_parameter_firstdataline)
elseif density_parameter_type == "BSFG"
    # Load DataFile
end

if isobaric_distribution_type == "MEAN_VALUES"
    dpAZ = DataFrame(A = NaN, rms_A = 0.6, ΔZ_A = -0.5)
elseif isobaric_distribution_type == "DATA"
    #Load dpAZ from file
end

#Revert relative PATH to project root folder
cd(@__DIR__)

println("*begin DSE computation")
#Begin program execution

fdmn = Fragmentation_domain(236, 92, 3, 118, 160, dpAZ)

println("*program execution successful")