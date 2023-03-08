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
include("dse_eq_solvers.jl")

println("*building fragmentation domain")
fragmdomain = Fragmentation_domain(A₀, Z₀, No_ZperA, A_H_min, A_H_max, isobaric_distribution_datafile)

println("*partitioning Total Excitation Energy")
E_excitation = TXE_partitioning(txe_partitioning_type, A₀, Z₀, A_H_min, A_H_max, Eₙ, fragmdomain, txe_partitioning_datafile, tkerange, density_parameter_type, density_parameter_datafile, dmass_excess)

println("*solving DSE energy conservation equations")
DSE_eq_output = DSE_equation_solver(evaporation_cs_type, E_excitation, density_parameter_type, density_parameter_datafile, dmass_excess)

println("*processing DSE equations primary output")
Raw_output_datafile = Process_main_output(DSE_eq_output, evaporation_cs_type)

if write_primary_output == "YES"
    println("*writing main DSE output data to file")
    Write_seq_output(A₀, Z₀, A_H_min, A_H_max, No_ZperA, Eₙ, tkerange, fragmdomain, E_excitation, Raw_output_datafile, density_parameter_type, density_parameter_datafile, evaporation_cs_type, mass_excess_filename, txe_partitioning_type, dmass_excess)
end

if secondary_outputs == "YES"
    println("*averaging data over $yield_distribution_filename experimental Yield distribution")
    include("secondary_outputs.jl")
end

if neutron_spectrum == "YES"
    println("*generating prompt neutron spectrum")
    include("neutron_spectrum.jl")
end

if generate_plots == "YES"
    println("*plotting data")
    include("plots.jl")
end

println("*program execution successful!")