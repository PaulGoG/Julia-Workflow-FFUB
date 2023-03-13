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
include("dse_eq_solvers.jl")

#Execute function calls
fragmdomain = Fragmentation_domain(A₀, Z₀, No_ZperA, A_H_range, isobaric_distribution_data)
E_excitation = TXE_partitioning(txe_partitioning_type, A₀, Z₀, A_H_range, Eₙ, fragmdomain, txe_partitioning_data, tkerange, density_parameter_type, density_parameter_data, dmass_excess)
DSE_eq_output = DSE_equation_solver(evaporation_cs_type, E_excitation, density_parameter_type, density_parameter_data, dmass_excess)
Raw_output_datafile = Process_main_output(DSE_eq_output, evaporation_cs_type)
if write_primary_output == "YES"
    Write_seq_output(A₀, Z₀, A_H_min, A_H_max, No_ZperA, Eₙ, tkerange, fragmdomain, E_excitation, Raw_output_datafile, density_parameter_type, density_parameter_data, fissionant_nucleus_identifier, mass_excess_filename, txe_partitioning_type, txe_partitioning_data, evaporation_cs_type, dmass_excess)
end
if secondary_outputs == "YES"
    include("secondary_outputs.jl")
end
if neutron_spectrum == "YES"
    include("neutron_spectrum.jl")
end
if generate_plots == "YES"
    include("plots.jl")
end
println("*end program execution at $(Dates.format(now(), "HH:MM:SS"))")
println("*program execution succesful!")