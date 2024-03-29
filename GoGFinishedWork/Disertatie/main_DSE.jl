#Main body of the DSE (Deterministic Sequential Emission) program
#####
include("init.jl")
include("input.jl")
include("options.jl")
include("error_handling.jl")
include("aux_func.jl")
include("density_parameters.jl")
include("txe_partitioning.jl")
include("dse_eq_solvers.jl")
include("write_output.jl")

fragmdomain = Fragmentation_domain(A₀, Z₀, No_ZperA, A_H_range, isobaric_distribution_data)
E_excitation = TXE_partitioning(txe_partitioning_type, A₀, Z₀, A_H_range, fission_type, E_incident, 
                fragmdomain, txe_partitioning_data, tkerange, density_parameter_type, density_parameter_data, dmass_excess)
DSE_eq_output = DSE_equation_solver(evaporation_cs_type, E_excitation, density_parameter_type, density_parameter_data, dmass_excess)
Raw_output_datafile = Process_main_output(DSE_eq_output, evaporation_cs_type)
if !secondary_outputs
    Write_seq_output_P(A₀, Z₀, A_H_min, A_H_max, No_ZperA, tkerange, 
    fragmdomain, E_excitation, Raw_output_datafile, density_parameter_type, density_parameter_data, 
    fissionant_nucleus_identifier, file_output_identifier, mass_excess_filename, 
    txe_partitioning_type, txe_partitioning_data, evaporation_cs_type, dmass_excess)
else 
    include("secondary_outputs.jl")
end
if neutron_spectrum
    include("neutron_spectrum.jl")
end
if generate_plots
    include("plots.jl")
end
println("*end program execution at $(Dates.format(now(), "HH:MM:SS"))")
println("*program execution succesful!")