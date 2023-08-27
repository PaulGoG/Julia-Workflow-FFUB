#Write output to file with p(Z,A) or Y(A,Z,TKE)
function Write_seq_output_P(A_0, Z_0, A_H_min, A_H_max, No_ZperA, tkerange, 
    fragmdomain, E_excitation, Processed_raw_output, density_parameter_type, density_parameter_data, 
    fissionant_nucleus_identifier, file_output_identifier, mass_excess_filename, 
    txe_partitioning_type, txe_partitioning_data, evaporation_cs_type, dm)
    
    println("*writing DSE output data to file")
    horizontal_delimiter = lpad('-', 159, '-')

    open("$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_$(file_output_identifier)_readme.out", "w") do file
        write(file, "DSE main output file generated at $(Dates.format(now(), "HH:MM:SS")) corresponding to input data:\r\n")
        write(file, "$(fissionant_nucleus_identifier) $(file_output_identifier) (A₀ = $A_0, Z₀ = $Z_0), fission type: $fission_type, $No_ZperA Z per A, mass excess file - $mass_excess_filename\r\n")
        write(file, "Heavy Fragment mass number ranges from $A_H_min to $A_H_max\r\n")
        write(file, "TKE ∈ $tkerange\r\n")
        write(file, "TXE partitioning method - $txe_partitioning_type\r\n")
        if txe_partitioning_type == "MSCZ"
            write(file, "TXE partitioning data used: Extra deformation energies from $txe_partitioning_filename\r\n")
        elseif txe_partitioning_type == "RT"
            write(file, "TXE partitioning data used: RT(A_H) denoted by segments $txe_partitioning_data\r\n")
        elseif txe_partitioning_type == "PARAM"
            write(file, "TXE partitioning data used: Ratio(A_H) = E*_H/TXE denoted by segments $txe_partitioning_data\r\n")
        end
        if evaporation_cs_type == "CONSTANT"
            write(file, "Neutron evaporation cross section is considered CONSTANT\r\n")
        elseif evaporation_cs_type == "VARIABLE"
            write(file, "Neutron evaporation cross section is considered VARIABLE and calculated using s-wave neutron force function\r\n")
        end
        write(file, "$horizontal_delimiter")
    end

    open("$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_$(file_output_identifier)_main_DSE.out", "w") do file
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                Z_L = Z_0 - Z_H
                P_Z_A = round(fragmdomain.Value[(fragmdomain.A .== A_H) .& (fragmdomain.Z .== Z_H)][1], digits = 8)
                a_L = round(density_parameter(density_parameter_type, A_L, Z_L, density_parameter_data), digits = 8)
                a_H = round(density_parameter(density_parameter_type, A_H, Z_H, density_parameter_data), digits = 8)
                for TKE in unique(Processed_raw_output.TKE[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H)])
                    E_excit_H = round(E_excitation.Value[(E_excitation.A .== A_H) .& (E_excitation.Z .== Z_H) .& (E_excitation.TKE .== TKE)][1], digits = 8)
                    E_excit_L = round(E_excitation.Value[(E_excitation.A .== A_L) .& (E_excitation.Z .== Z_L) .& (E_excitation.TKE .== TKE)][1], digits = 8)
                    write(file, "$A_H $Z_H $TKE $P_Z_A $E_excit_L $E_excit_H $a_L $a_H\r\n")
                        
                    n_range_L = Processed_raw_output.No_Sequence[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE)]
                    if last(n_range_L) != 0
                        write(file, "$(last(n_range_L))\r\n")
                        for k in n_range_L
                            S_k = Separation_energy(1, 0, A_L - k + 1, Z_L, dm)[1]
                            write(file, "$(round(S_k, digits = 8)) ")
                        end
                        write(file, "\r\n")
                        for k in n_range_L
                            a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            write(file, "$(round(a_k, digits = 8)) ")
                        end
                        write(file, "\r\n")
                        for k in n_range_L
                            T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            if T_k > 0 
                                write(file, "$(round(T_k, digits = 8)) ")
                            else
                                write(file, "0.00000000 ")
                            end
                        end
                        write(file, "\r\n")
                        for k in n_range_L
                            T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            if T_k > 0 
                                a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                Eᵣ_k = Energy_FermiGas(a_k, T_k)
                                write(file, "$(round(Eᵣ_k, digits = 8)) ")
                            else
                                write(file, "0.00000000 ")
                            end
                        end
                        write(file, "\r\n")
                        for k in n_range_L
                            avgε_k = Processed_raw_output.Avg_εₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            if avgε_k > 0 
                                write(file, "$(round(avgε_k, digits = 8)) ")
                            else
                                write(file, "0.00000000 ")
                            end
                        end
                        write(file, "\r\n")
                    else
                        write(file, "0\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n")
                    end

                    n_range_H = Processed_raw_output.No_Sequence[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE)]
                    if last(n_range_H) != 0
                        write(file, "$(last(n_range_H))\r\n")
                        for k in n_range_H
                            S_k = Separation_energy(1, 0, A_H - k + 1, Z_H, dm)[1]
                            write(file, "$(round(S_k, digits = 8)) ")
                        end
                        write(file, "\r\n")
                        for k in n_range_H
                            a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            write(file, "$(round(a_k, digits = 8)) ")
                        end
                        write(file, "\r\n")
                        for k in n_range_H
                            T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            if T_k > 0 
                                write(file, "$(round(T_k, digits = 8)) ")
                            else
                                write(file, "0.00000000 ")
                            end
                        end
                        write(file, "\r\n")
                        for k in n_range_H
                            T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            if T_k > 0 
                                a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                Eᵣ_k = Energy_FermiGas(a_k, T_k)
                                write(file, "$(round(Eᵣ_k, digits = 8)) ")
                            else
                                write(file, "0.00000000 ")
                            end
                        end
                        write(file, "\r\n")
                        for k in n_range_H
                            avgε_k = Processed_raw_output.Avg_εₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            if avgε_k > 0 
                                write(file, "$(round(avgε_k, digits = 8)) ")
                            else
                                write(file, "0.00000000 ")
                            end
                        end
                        write(file, "\r\n")
                    else
                        write(file, "0\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n")
                    end
                end
            end
        end
    end
end
function Write_seq_output_Y(A_0, Z_0, A_H_min, A_H_max, No_ZperA, tkerange, 
    y_A_Z_TKE, E_excitation, Processed_raw_output, density_parameter_type, density_parameter_data, 
    fissionant_nucleus_identifier, file_output_identifier, mass_excess_filename, 
    txe_partitioning_type, txe_partitioning_data, evaporation_cs_type, dm)

    println("*writing DSE output data to file")
    horizontal_delimiter = lpad('-', 159, '-')

    open("$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_$(file_output_identifier)_readme.out", "w") do file
        write(file, "DSE main output file generated at $(Dates.format(now(), "HH:MM:SS")) corresponding to input data:\r\n")
        write(file, "$(fissionant_nucleus_identifier) $(file_output_identifier) (A₀ = $A_0, Z₀ = $Z_0), fission type: $fission_type, $No_ZperA Z per A, mass excess file - $mass_excess_filename\r\n")
        write(file, "Heavy Fragment mass number ranges from $A_H_min to $A_H_max\r\n")
        write(file, "TKE ∈ $tkerange\r\n")
        write(file, "TXE partitioning method - $txe_partitioning_type\r\n")
        if txe_partitioning_type == "MSCZ"
            write(file, "TXE partitioning data used: Extra deformation energies from $txe_partitioning_filename\r\n")
        elseif txe_partitioning_type == "RT"
            write(file, "TXE partitioning data used: RT(A_H) denoted by segments $txe_partitioning_data\r\n")
        elseif txe_partitioning_type == "PARAM"
            write(file, "TXE partitioning data used: Ratio(A_H) = E*_H/TXE denoted by segments $txe_partitioning_data\r\n")
        end
        if evaporation_cs_type == "CONSTANT"
            write(file, "Neutron evaporation cross section is considered CONSTANT\r\n")
        elseif evaporation_cs_type == "VARIABLE"
            write(file, "Neutron evaporation cross section is considered VARIABLE and calculated using s-wave neutron force function\r\n")
        end
        write(file, "$horizontal_delimiter")
    end

    open("$(file_output_identifier)_output_data/$(fissionant_nucleus_identifier)_$(file_output_identifier)_main_DSE.out", "w") do file
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for Z_H in sort(unique(y_A_Z_TKE.Z[y_A_Z_TKE.A .== A_H]))
                Z_L = Z_0 - Z_H
                a_L = round(density_parameter(density_parameter_type, A_L, Z_L, density_parameter_data), digits = 8)
                a_H = round(density_parameter(density_parameter_type, A_H, Z_H, density_parameter_data), digits = 8)
                for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H)]
                    Y_A_Z_TKE = round(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H) .& (y_A_Z_TKE.TKE .== TKE)][1], digits = 8)
                    if isassigned(E_excitation.Value[(E_excitation.A .== A_H) .& (E_excitation.Z .== Z_H) .& (E_excitation.TKE .== TKE)], 1)
                        E_excit_H = round(E_excitation.Value[(E_excitation.A .== A_H) .& (E_excitation.Z .== Z_H) .& (E_excitation.TKE .== TKE)][1], digits = 8)
                        E_excit_L = round(E_excitation.Value[(E_excitation.A .== A_L) .& (E_excitation.Z .== Z_L) .& (E_excitation.TKE .== TKE)][1], digits = 8)
                        write(file, "$A_H $Z_H $TKE $Y_A_Z_TKE $E_excit_L $E_excit_H $a_L $a_H\r\n")
                        
                        n_range_L = Processed_raw_output.No_Sequence[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE)]
                        if last(n_range_L) != 0
                            write(file, "$(last(n_range_L))\r\n")
                            for k in n_range_L
                                S_k = Separation_energy(1, 0, A_L - k + 1, Z_L, dm)[1]
                                write(file, "$(round(S_k, digits = 8)) ")
                            end
                            write(file, "\r\n")
                            for k in n_range_L
                                a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                write(file, "$(round(a_k, digits = 8)) ")
                            end
                            write(file, "\r\n")
                            for k in n_range_L
                                T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                if T_k > 0 
                                    write(file, "$(round(T_k, digits = 8)) ")
                                else
                                    write(file, "0.00000000 ")
                                end
                            end
                            write(file, "\r\n")
                            for k in n_range_L
                                T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                if T_k > 0 
                                    a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                    Eᵣ_k = Energy_FermiGas(a_k, T_k)
                                    write(file, "$(round(Eᵣ_k, digits = 8)) ")
                                else
                                    write(file, "0.00000000 ")
                                end
                            end
                            write(file, "\r\n")
                            for k in n_range_L
                                avgε_k = Processed_raw_output.Avg_εₖ[(Processed_raw_output.A .== A_L) .& (Processed_raw_output.Z .== Z_L) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                if avgε_k > 0 
                                    write(file, "$(round(avgε_k, digits = 8)) ")
                                else
                                    write(file, "0.00000000 ")
                                end
                            end
                            write(file, "\r\n")
                        else
                            write(file, "0\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n")
                        end

                        n_range_H = Processed_raw_output.No_Sequence[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE)]
                        if last(n_range_H) != 0
                            write(file, "$(last(n_range_H))\r\n")
                            for k in n_range_H
                                S_k = Separation_energy(1, 0, A_H - k + 1, Z_H, dm)[1]
                                write(file, "$(round(S_k, digits = 8)) ")
                            end
                            write(file, "\r\n")
                            for k in n_range_H
                                a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                write(file, "$(round(a_k, digits = 8)) ")
                            end
                            write(file, "\r\n")
                            for k in n_range_H
                                T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                if T_k > 0 
                                    write(file, "$(round(T_k, digits = 8)) ")
                                else
                                    write(file, "0.00000000 ")
                                end
                            end
                            write(file, "\r\n")
                            for k in n_range_H
                                T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                if T_k > 0 
                                    a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                    Eᵣ_k = Energy_FermiGas(a_k, T_k)
                                    write(file, "$(round(Eᵣ_k, digits = 8)) ")
                                else
                                    write(file, "0.00000000 ")
                                end
                            end
                            write(file, "\r\n")
                            for k in n_range_H
                                avgε_k = Processed_raw_output.Avg_εₖ[(Processed_raw_output.A .== A_H) .& (Processed_raw_output.Z .== Z_H) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                                if avgε_k > 0 
                                    write(file, "$(round(avgε_k, digits = 8)) ")
                                else
                                    write(file, "0.00000000 ")
                                end
                            end
                            write(file, "\r\n")
                        else
                            write(file, "0\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n0.00000000\r\n")
                        end
                    end
                end
            end
        end
    end
end