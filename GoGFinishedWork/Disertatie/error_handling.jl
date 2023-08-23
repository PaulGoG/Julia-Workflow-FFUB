#Flow control and error handling for input parameters and filenames of the main program
#####
if A_H_min < A₀/2 || !isinteger(A_H_min)
    error("A_H_min must be an Integer greater than $(Int(round(A₀/2)))!")
end

if A_H_max >= A₀ || !isinteger(A_H_max)
    error("A_H_max must be an Integer lesser than $(A₀)!")
end

if A_H_min >= A_H_max
    error("invalid Heavy Fragment range!")
end

if TKE_min >= TKE_max 
    error("invalid TKE range!")
end

if !isfile(mass_excess_filename)
    error("$mass_excess_filename does not exist at input_data/  PATH!")
end

if isassigned(filter(x -> !in(x, ("Z", "A", "Symbol", "D", "σ_D")), mass_excess_header), 1)
    error(
    "mass_excess_header contains invalid fields: $(filter(x -> !in(x, ("Z", "A", "Symbol", "D", "σ_D")) ,mass_excess_header));
    allowed header names list: Z, A, Symbol, D, σ_D"
    )
end

if fission_type != "SF" && fission_type != "(n,f)"
    error("$fission_type is not a valid input for fission_type!")
end

if density_parameter_type != "GC" && density_parameter_type != "BSFG"
    error("$density_parameter_type is not a valid input for density_parameter_type!")
end

if density_parameter_type == "GC"
    if !isfile(density_parameter_filename)
        error("$density_parameter_filename does not exist at input_data/  PATH!")
    end
    if isassigned(filter(x -> !in(x, ("n", "S_Z", "S_N")), density_parameter_header), 1)
        error(
        "density_parameter_header contains invalid fields: $(filter(x -> !in(x, ("n", "S_Z", "S_N")), density_parameter_header));
        allowed header names list: n, S_Z, S_N"
        )
    end
end

if evaporation_cs_type != "CONSTANT" && evaporation_cs_type != "VARIABLE"
    error("$evaporation_cs_type is not a valid input for evaporation_cs_type!")
end

if isobaric_distribution_type != "MEAN_VALUES" && isobaric_distribution_type != "DATA"
    error("$isobaric_distribution_type is not a valid input for isobaric_distribution_type!")
end

if isobaric_distribution_type == "DATA"
    if !isfile(isobaric_distribution_filename)
        error("$isobaric_distribution_filename does not exist at input_data/  PATH!")
    end
    if isassigned(filter(x -> !in(x, ("A", "ΔZ_A", "rms_A")), isobaric_distribution_header), 1)
        error(
        "isobaric_distribution_header contains invalid fields: $(filter(x -> !in(x, ("A", "ΔZ_A", "rms_A")), isobaric_distribution_header));
        allowed header names list: A, ΔZ_A, rms_A"
        )
    end
end

if !isodd(No_ZperA)
    error("No_ZperA must be an odd Integer!")
end

if txe_partitioning_type != "MSCZ" && txe_partitioning_type != "PARAM" && txe_partitioning_type != "RT"
    error("$txe_partitioning_type is not a valid input for txe_partitioning_type!")
end

if txe_partitioning_type == "MSCZ"
    if !isfile(txe_partitioning_filename)
        error("$txe_partitioning_filename does not exist at input_data/  PATH!")
    end
    if isassigned(filter(x -> !in(x, ("A", "Z", "Value")), txe_partitioning_header), 1)
        error(
        "txe_partitioning_header contains invalid fields: $(filter(x -> !in(x, ("A", "Z", "Value")), txe_partitioning_header));
        allowed header names list: A, Z, Value"
        )
    end
end

if !secondary_outputs
    if secondary_output_Yield
        error("averaged yield distributions cannot be calculated without valid yield distribution!")
    elseif secondary_output_TXE_Q
        error("Q-value and TXE averaged distributions cannot be calculated without valid yield distribution!")
    elseif secondary_output_E_excitation
        error("E* averaged distributions cannot be calculated without valid yield distribution!")
    elseif secondary_output_nu
        error("averaged neutron multiplicities cannot be calculated without valid yield distribution!")
    elseif secondary_output_Ap
        error("post emission distributions cannot be calculated without valid yield distribution!")
    elseif secondary_output_T
        error("averaged residual temperatures cannot be calculated without valid yield distribution!")
    elseif secondary_output_Tₖ
        error("averaged residual temperatures cannot be calculated without valid yield distribution!")
    elseif secondary_output_avg_ε
        error("averaged center-of-mass neutron energies cannot be calculated without valid yield distribution!")
    elseif secondary_output_avg_εₖ
        error("averaged center-of-mass neutron energies cannot be calculated without valid yield distribution!")
    elseif secondary_output_Eᵣ
        error("averaged residual fragment energies cannot be calculated without valid yield distribution!")
    elseif generate_plots
        error("plots cannot be generated without yield-averaged quantities!")
    elseif neutron_spectrum
        error("prompt neutron spectrum cannot be calculated without valid yield distribution!")
    end
else
    if !isfile(yield_distribution_filename)
        error("$yield_distribution_filename does not exist at input_data/  PATH!")
    end
    if isassigned(filter(x -> !in(x, ("A", "TKE", "Value", "σ")), yield_distribution_header), 1)
        error(
        "yield_distribution_header contains invalid fields: $(filter(x -> !in(x, ("A", "TKE", "Value", "σ")), yield_distribution_header));
        allowed header names list: A, TKE, Value, σ"
        )
    end
    if neutron_spectrum && E_min >= E_max
        error("invalid neutron spectrum Energy range!")
    end
end