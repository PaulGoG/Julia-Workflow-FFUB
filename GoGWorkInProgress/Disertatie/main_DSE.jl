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

#Revert relative PATH to project root folder
cd(@__DIR__)

println("*begin DSE computation")
println("*building fragmentation domain")
fragmdomain = Fragmentation_domain(A₀, Z₀, No_ZperA, A_H_min, A_H_max, isobaric_distribution_datafile)

println("*partitioning Total Excitation Energy")
E_excitation = TXE_partitioning(txe_partitioning_type, A₀, Z₀, A_H_min, A_H_max, Eₙ, fragmdomain, txe_partitioning_datafile, tkerange, density_parameter_type, density_parameter_datafile, dmass_excess)

println("*solving DSE energy conservation equations")
DSE_eq_output = DSE_equation_solver(evaporation_cs_type, E_excitation, density_parameter_type, density_parameter_datafile, dmass_excess)

println("*preparing output datafile")
Raw_output_datafile = Process_main_output(DSE_eq_output, evaporation_cs_type)

println("*writing output to file")
Write_seq_output(A₀, Z₀, No_ZperA, Eₙ, E_excitation, Raw_output_datafile, density_parameter_type, density_parameter_datafile, evaporation_cs_type, dmass_excess)

println("*main program execution successful!")

#=
ν_A_Z_TKE = Neutron_multiplicity_A_Z_TKE(DataFrame(
    A = Raw_output_datafile.A,
    Z = Raw_output_datafile.Z,
    TKE = Raw_output_datafile.TKE,
    No_Sequence = Raw_output_datafile.No_Sequence
))

T_A_Z_TKE = SeqAvg_A_Z_TKE(DataFrame(
    A = Raw_output_datafile.A,
    Z = Raw_output_datafile.Z,
    TKE = Raw_output_datafile.TKE,
    No_Sequence = Raw_output_datafile.No_Sequence,
    Value = Raw_output_datafile.Tₖ
))
=#
