function E_f_L(A_0, A_L, A_H, TKE)
    return (A_H/A_L)*TKE/A_0
end
function E_f_H(A_0, A_L, A_H, TKE)
    return (A_L/A_H)*TKE/A_0
end
function u_1(E, E_f)
    return (sqrt(E) - sqrt(E_f))^2 
end
function u_2(E, E_f)
    return (sqrt(E) + sqrt(E_f))^2
end
function Neutron_spectrum_CONSTANT_cs(A_0, Z_0, A_H_min, A_H_max, energyrange, Raw_output_datafile::DataFrame, y_A_Z_TKE::DataFrame)
    n_E = Neutron_spectrum_object(Float64[], Float64[])
    function Neutron_spectrum_k(ε::Float64, T::Float64)
        return ε *exp(-ε/T) /T^2
    end
    E_min = first(energyrange)
    ΔE = last(energyrange) - first(energyrange)
    for E in energyrange
        Denominator = 0.0
        Numerator = 0.0
        if isinteger(E - E_min)
            println("Progress: $(Int(round(100*(E - E_min)/ΔE)))%")
        end
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for Z_H in unique(y_A_Z_TKE.Z[(y_A_Z_TKE.A .== A_H)])
                Z_L = Z_0 - Z_H
                for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H)]
                    Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H) .& (y_A_Z_TKE.TKE .== TKE)][1]
                    if A_L == A_H
                        Y_A_Z_TKE /= 2
                    end
                    #HF
                    if !isassigned(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)], 1)
                        ν_H = 0.0
                        N_H = 0.0
                    elseif last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]) == 0
                        ν_H = 0.0
                        N_H = 0.0                    
                    else 
                        ν_H = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)])
                        T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]
                        E_f = E_f_H(A_0, A_L, A_H, TKE)
                        Integrand_H(ε::Float64) =(sum(Neutron_spectrum_k.(ε, T)) /ν_H)/(4 *sqrt(ε*E_f))
                        N_H = quadgk(Integrand_H, u_1(E, E_f), u_2(E, E_f))[1]
                    end
                    #LF
                    if !isassigned(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)], 1)
                        ν_L = 0.0
                        N_L = 0.0
                    elseif last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]) == 0
                        ν_L = 0.0
                        N_L = 0.0
                    else 
                        ν_L = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)])
                        T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]
                        E_f = E_f_L(A_0, A_L, A_H, TKE)
                        Integrand_L(ε::Float64) = (sum(Neutron_spectrum_k.(ε, T)) /ν_L)/(4 *sqrt(ε*E_f))
                        N_L = quadgk(Integrand_L, u_1(E, E_f), u_2(E, E_f))[1]
                    end
                    #Add normalised result contribution to spectrum
                    ν_total = ν_H + ν_L
                    if ν_total > 0
                        Denominator += Y_A_Z_TKE
                        Numerator += (N_H *ν_H + N_L *ν_L) *Y_A_Z_TKE /ν_total
                    end
                end
            end
        end
        if Denominator > 0
            push!(n_E.E, E)
            push!(n_E.Value, Numerator/Denominator)
        end
    end
    return n_E
end
function Neutron_spectrum_VARIABLE_cs(A_0, Z_0, A_H_min, A_H_max, energyrange, Raw_output_datafile::DataFrame, y_A_Z_TKE::DataFrame)
    n_E = Neutron_spectrum_object(Float64[], Float64[])
    function Neutron_spectrum_k(ε::Float64, α::Float64, T::Float64)
        (ε + α*sqrt(ε)) *exp(-ε/T) /((sqrt(T) +α*sqrt(π)/2) *T^(3/2))
    end
    E_min = first(energyrange)
    ΔE = last(energyrange) - first(energyrange)
    for E in energyrange
        Denominator = 0.0
        Numerator = 0.0
        if isinteger(E - E_min)
            println("Progress: $(Int(round(100*(E - E_min)/ΔE)))%")
        end
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for Z_H in unique(y_A_Z_TKE.Z[(y_A_Z_TKE.A .== A_H)])
                Z_L = Z_0 - Z_H
                for TKE in y_A_Z_TKE.TKE[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H)]
                    Y_A_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H) .& (y_A_Z_TKE.TKE .== TKE)][1]
                    if A_L == A_H
                        Y_A_Z_TKE /= 2
                    end
                    #HF
                    if !isassigned(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)], 1)
                        ν_H = 0.0
                        N_H = 0.0
                    elseif last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]) == 0
                        ν_H = 0.0
                        N_H = 0.0                    
                    else 
                        ν_H = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)])
                        T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]
                        α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]
                        E_f = E_f_H(A_0, A_L, A_H, TKE)
                        Integrand_H(ε::Float64) =(sum(Neutron_spectrum_k.(ε, α, T)) /ν_H)/(4 *sqrt(ε*E_f))
                        N_H = quadgk(Integrand_H, u_1(E, E_f), u_2(E, E_f))[1]
                    end
                    #LF
                    if !isassigned(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)], 1)
                        ν_L = 0.0
                        N_L = 0.0
                    elseif last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]) == 0
                        ν_L = 0.0
                        N_L = 0.0
                    else 
                        ν_L = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)])
                        T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]
                        α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]
                        E_f = E_f_L(A_0, A_L, A_H, TKE)
                        Integrand_L(ε::Float64) = (sum(Neutron_spectrum_k.(ε, α, T)) /ν_L)/(4 *sqrt(ε*E_f))
                        N_L = quadgk(Integrand_L, u_1(E, E_f), u_2(E, E_f))[1]
                    end
                    #Add normalised result contribution to spectrum
                    ν_total = ν_H + ν_L
                    if ν_total > 0
                        Denominator += Y_A_Z_TKE
                        Numerator += (N_H *ν_H + N_L *ν_L) *Y_A_Z_TKE /ν_total
                    end
                end
            end
        end
        if Denominator > 0
            push!(n_E.E, E)
            push!(n_E.Value, Numerator/Denominator)
        end
    end
    return n_E
end
function Neutron_spectrum_builder(A_0, Z_0, A_H_min, A_H_max, energyrange, Raw_output_datafile::DataFrame, evaporation_cs_type, y_A_Z_TKE::DataFrame)
    println("*begin generating prompt neutron spectrum at $(Dates.format(now(), "HH:MM:SS"))")
    if evaporation_cs_type == "CONSTANT"
        return Neutron_spectrum_CONSTANT_cs(A_0, Z_0, A_H_min, A_H_max, energyrange, Raw_output_datafile, y_A_Z_TKE)
    elseif evaporation_cs_type == "VARIABLE"
        return Neutron_spectrum_VARIABLE_cs(A_0, Z_0, A_H_min, A_H_max, energyrange, Raw_output_datafile, y_A_Z_TKE)
    end
end
function Normalise_spectrum_to_Maxwellian(E::Vector{Float64}, N::Vector{Float64}, T_M::Float64)
    Maxwell(ε::Float64) = (2/sqrt(π)) *T_M^(-3/2) *sqrt(ε) *exp(-ε/T_M)
    Surface_Maxwell = quadgk(Maxwell, first(E), last(E))[1]
    Surface_spectrum = trapz(E, N)
    norm_parameter = Surface_Maxwell/Surface_spectrum
    for index in eachindex(E)
       N[index] *= norm_parameter/Maxwell(E[index])
    end
end
#####

Adjusted_yield = DataFrame(
    A = y_A_Z_TKE.A[(y_A_Z_TKE.Value .>= Yield_cutoff_value) .& (y_A_Z_TKE.A .>= A_H_min)], 
    Z = y_A_Z_TKE.Z[(y_A_Z_TKE.Value .>= Yield_cutoff_value) .& (y_A_Z_TKE.A .>= A_H_min)],
    TKE = y_A_Z_TKE.TKE[(y_A_Z_TKE.Value .>= Yield_cutoff_value) .& (y_A_Z_TKE.A .>= A_H_min)],
    Value = y_A_Z_TKE.Value[(y_A_Z_TKE.Value .>= Yield_cutoff_value) .& (y_A_Z_TKE.A .>= A_H_min)]
)
n_E = Neutron_spectrum_builder(A₀, Z₀, A_H_min, A_H_max, energyrange, Raw_output_datafile, evaporation_cs_type, Adjusted_yield)
T_M_eq = 2/3 *trapz(n_E.E, n_E.Value .* n_E.E)/trapz(n_E.E, n_E.Value)
Ratio_to_Maxwellian = copy(n_E.Value)
Normalise_spectrum_to_Maxwellian(n_E.E, Ratio_to_Maxwellian, T_M_eq)
CSV.write(
    "output_data/$(fissionant_nucleus_identifier)_neutron_spectrum_TMequiv_$(round(T_M_eq, digits = 4)).OUT", 
    DataFrame(E = n_E.E, N = n_E.Value, Ratio_Maxwellian = Ratio_to_Maxwellian), 
    writeheader=true, newline="\r\n", delim=' '
)