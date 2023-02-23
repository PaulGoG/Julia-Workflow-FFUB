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
    density_parameter_datafile = CSV.read(density_parameter_filename, DataFrame; delim = density_parameter_delimiter, ignorerepeated = true, header = density_parameter_header, skipto = density_parameter_firstdataline)
elseif density_parameter_type == "BSFG"
    density_parameter_datafile = dm
end

if isobaric_distribution_type == "MEAN_VALUES"
    dpAZ = DataFrame(A = NaN, rms_A = NaN, ΔZ_A = NaN)
elseif isobaric_distribution_type == "DATA"
    #Load dpAZ from file
end

if txe_partitioning_type == "MSCZ"
    txe_partitioning_datafile = CSV.read(txe_partitioning_filename, DataFrame; delim = txe_partitioning_delimiter, ignorerepeated = true, header = txe_partitioning_header, skipto = txe_partitioning_firstdataline)
elseif txe_partitioning_type == "PARAM"
    txe_partitioning_datafile = CSV.read(txe_partitioning_filename, DataFrame; delim = txe_partitioning_delimiter, ignorerepeated = true, header = txe_partitioning_header, skipto = txe_partitioning_firstdataline)
end

#Revert relative PATH to project root folder
cd(@__DIR__)

println("*begin DSE computation")
#Begin program execution

println("*building fragmentation domain")
fragmdomain = Fragmentation_domain(A₀, Z₀, No_ZperA, A_H_min, A_H_max, dpAZ)

println("*partitioning Total Excitation Energy")
E_excitation = TXE_partitioning(txe_partitioning_type, A₀, Z₀, A_H_min, A_H_max, Eₙ, fragmdomain, txe_partitioning_datafile, tkerange, density_parameter_type, density_parameter_datafile, dm)

println(E_excitation)

println("*solving DSE equations")

println("*preparing output datafile")

println("*program execution successful")