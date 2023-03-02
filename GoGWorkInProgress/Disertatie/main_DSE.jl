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
DSE_eq_output = DSE_equation_solver(evaporation_cs_type, fragmdomain, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dmass_excess)

println("*preparing output datafile")
Raw_output_datafile = Process_main_output(DSE_eq_output, evaporation_cs_type)
#Write_seq_output(Raw_output_datafile, evaporation_cs_type)

println("*main program execution successful!")

#=
Fisiere separate de output post-mediere pe secvente!

ν_A_Z_TKE = Neutron_multiplicity_A_Z_TKE(Output_datafile);

T_A_Z_TKE = SeqAvg_A_Z_TKE(DataFrame(
    A = Output_datafile.A,
    Z = Output_datafile.Z,
    TKE = Output_datafile.TKE,
    No_Sequence = Output_datafile.No_Sequence,
    Value = Output_datafile.Tₖ
));
=#
