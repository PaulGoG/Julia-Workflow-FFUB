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
include("secondary_outputs.jl")

#Revert relative PATH to project root folder
cd(@__DIR__)

println("*begin DSE computation")
println("*building fragmentation domain")
fragmdomain = Fragmentation_domain(A₀, Z₀, No_ZperA, A_H_min, A_H_max, isobaric_distribution_datafile)

println("*partitioning Total Excitation Energy")
E_excitation = TXE_partitioning(txe_partitioning_type, A₀, Z₀, A_H_min, A_H_max, Eₙ, fragmdomain, txe_partitioning_datafile, tkerange, density_parameter_type, density_parameter_datafile, dmass_excess)

println("*solving DSE energy conservation equations")
DSE_eq_output = DSE_equation_solver(evaporation_cs_type, E_excitation, density_parameter_type, density_parameter_datafile, dmass_excess)

println("*processing equation output")
Raw_output_datafile = Process_main_output(DSE_eq_output, evaporation_cs_type)

if write_primary_outputs == "YES"
    println("*writing output to file")
    Write_seq_output(A₀, Z₀, A_H_min, A_H_max, No_ZperA, Eₙ, tkerange, fragmdomain, E_excitation, Raw_output_datafile, density_parameter_type, density_parameter_datafile, evaporation_cs_type, mass_excess_filename, txe_partitioning_type, dmass_excess)
end

println("*main program execution successful!")

if secondary_outputs == "YES"
    println("*averaging data over $yield_distribution_filename experimental Yield distribution")
    y_A_Z_TKE = Process_yield_data(A₀, fragmdomain, dY)
    ν_A_Z_TKE = Neutron_multiplicity_A_Z_TKE(DataFrame(
    A = Raw_output_datafile.A,
    Z = Raw_output_datafile.Z,
    TKE = Raw_output_datafile.TKE,
    No_Sequence = Raw_output_datafile.No_Sequence
    ))
    ν_A = Average_over_TKE_Z(ν_A_Z_TKE, y_A_Z_TKE)
    Output = DataFrame(A = ν_A.Argument, ν = ν_A.Value)
    CSV.write("output_data/nu_A.OUT", Output, writeheader=true, newline='\n', delim=' ')
    ν_TKE = Average_over_A_Z(ν_A_Z_TKE, y_A_Z_TKE)
    Output = DataFrame(TKE = ν_TKE.Argument, ν = ν_TKE.Value)
    CSV.write("output_data/nu_TKE.OUT", Output, writeheader=true, newline='\n', delim=' ')
    y_Ap = Yield_post_neutron(y_A_Z_TKE, ν_A)
    Output = DataFrame(Aₚ = y_Ap.Argument, Y = y_Ap.Value, σ = y_Ap.σ)
    CSV.write("output_data/Y_Ap.OUT", Output, writeheader=true, newline='\n', delim=' ')
    y_Z_Ap = Yield_post_neutron(y_A_Z_TKE, ν_A_Z_TKE)
    Output = DataFrame(Aₚ = y_Z_Ap.A, Z = y_Z_Ap.Z, Y = y_Z_Ap.Value, σ = y_Z_Ap.σ)
    CSV.write("output_data/Y_Z_Ap.OUT", Output, writeheader=true, newline='\n', delim=' ')
end

#Testing for plots
Grid = zeros(length(unique(y_Z_Ap.A)), length(unique(y_Z_Ap.Z)))
i1 = 0
for A in unique(y_Z_Ap.A)
    i1 += 1
    i2 = 1
    for Z in unique(y_Z_Ap.Z)
        if isassigned(y_Z_Ap.Value[(y_Z_Ap.A .== A) .& (y_Z_Ap.Z .== Z)], 1)
            Grid[i1, i2] = y_Z_Ap.Value[(y_Z_Ap.A .== A) .& (y_Z_Ap.Z .== Z)][1]  
        end
        i2 += 1
    end
end
using Plots
surface(Grid')

#=
T_A_Z_TKE = SeqAvg_A_Z_TKE(DataFrame(
    A = Raw_output_datafile.A,
    Z = Raw_output_datafile.Z,
    TKE = Raw_output_datafile.TKE,
    No_Sequence = Raw_output_datafile.No_Sequence,
    Value = Raw_output_datafile.Tₖ
))
=#
#if plots

#if Neutron_spectrum