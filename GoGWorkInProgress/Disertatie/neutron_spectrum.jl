#Iteratively compute total neutron spectrum in Laboratory frame from main DSE data
#####
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
function Neutron_spectrum_CONSTANT_cs(A_0, Z_0, A_H_min, A_H_max, fragmdomain::Distribution, tkerange, energyrange, Raw_output_datafile::DataFrame, y_A_Z_TKE::Distribution)
    n_E = Neutron_spectrum_obj(Float64[], Float64[])
    function Neutron_spectrum_k(ε::Float64, T::Float64)
        return ε*exp(-ε/T)/T^2
    end
    for E in energyrange
        Denominator = 0.0
        Numerator = 0.0
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                Z_L = Z_0 - Z_H
                for TKE in tkerange
                    #HF
                    if !isassigned(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)], 1)
                        w_H = 0.0
                        N_H = 0.0
                        Y_AH_Z_TKE = 0.0
                    elseif last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]) == 0
                        w_H = 0.0
                        N_H = 0.0
                        Y_AH_Z_TKE = 0.0                       
                    elseif !isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        w_H = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)])
                        N_H = 0.0
                        Y_AH_Z_TKE = 0.0
                    else 
                        Y_AH_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        w_H = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)])
                        T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]
                        E_f = E_f_H(A_0, A_L, A_H, TKE)
                        u_1_H = u_1(E, E_f)
                        u_2_H = u_2(E, E_f)
                        Integrand_H(ε::Float64) =(sum(Neutron_spectrum_k.(ε, α, T)) /w_H)/(4 *sqrt(ε*E_f))
                        N_H = quadgk(Integrand_H, u_1_H, u_2_H)[1]
                    end
                    #LF
                    if !isassigned(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)], 1)
                        w_L = 0.0
                        N_L = 0.0
                        Y_AL_Z_TKE = 0.0
                    elseif last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]) == 0
                        w_L = 0.0
                        N_L = 0.0
                        Y_AL_Z_TKE = 0.0                       
                    elseif !isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_L) .& (y_A_Z_TKE.Z .== Z_L) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        w_L = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)])
                        N_L = 0.0
                        Y_AL_Z_TKE = 0.0
                    else 
                        Y_AL_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_L) .& (y_A_Z_TKE.Z .== Z_L) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        w_L = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)])
                        T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]
                        E_f = E_f_L(A_0, A_L, A_H, TKE)
                        u_1_L = u_1(E, E_f)
                        u_2_L = u_2(E, E_f)
                        Integrand_L(ε::Float64) = (sum(Neutron_spectrum_k.(ε, T)) /w_L)/(4 *sqrt(ε*E_f))
                        N_L = quadgk(Integrand_L, u_1_L, u_2_L)[1]
                    end
                    Denominator += Y_AH_Z_TKE + Y_AL_Z_TKE
                    Numerator += Y_AH_Z_TKE*w_H*N_H + Y_AL_Z_TKE*w_L*N_L
                end
            end
        end
        if Denominator != 0
            push!(n_E.E, E)
            push!(n_E.Value, Numerator/Denominator)
        end
    end
    return n_E
end
function Neutron_spectrum_VARIABLE_cs(A_0, Z_0, A_H_min, A_H_max, fragmdomain::Distribution, tkerange, energyrange, Raw_output_datafile::DataFrame, y_A_Z_TKE::Distribution)
    n_E = Neutron_spectrum_obj(Float64[], Float64[])
    function Neutron_spectrum_k(ε::Float64, α::Float64, T::Float64)
        (ε + α*sqrt(ε))*exp(-ε/T) /((sqrt(T) +α*sqrt(π)/2) *T^(3/2))
    end
    for E in energyrange
        Denominator = 0.0
        Numerator = 0.0
        for A_H in A_H_min:A_H_max
            A_L = A_0 - A_H
            for Z_H in fragmdomain.Z[fragmdomain.A .== A_H]
                Z_L = Z_0 - Z_H
                for TKE in tkerange
                    #HF
                    if !isassigned(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)], 1)
                        w_H = 0.0
                        N_H = 0.0
                        Y_AH_Z_TKE = 0.0
                    elseif last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]) == 0
                        w_H = 0.0
                        N_H = 0.0
                        Y_AH_Z_TKE = 0.0                       
                    elseif !isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        w_H = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)])
                        N_H = 0.0
                        Y_AH_Z_TKE = 0.0
                    else 
                        Y_AH_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_H) .& (y_A_Z_TKE.Z .== Z_H) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        w_H = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)])
                        T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]
                        α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]
                        E_f = E_f_H(A_0, A_L, A_H, TKE)
                        u_1_H = u_1(E, E_f)
                        u_2_H = u_2(E, E_f)
                        Integrand_H(ε::Float64) =(sum(Neutron_spectrum_k.(ε, α, T)) /w_H)/(4 *sqrt(ε*E_f))
                        N_H = quadgk(Integrand_H, u_1_H, u_2_H)[1]
                    end
                    #LF
                    if !isassigned(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)], 1)
                        w_L = 0.0
                        N_L = 0.0
                        Y_AL_Z_TKE = 0.0
                    elseif last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]) == 0
                        w_L = 0.0
                        N_L = 0.0
                        Y_AL_Z_TKE = 0.0                       
                    elseif !isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_L) .& (y_A_Z_TKE.Z .== Z_L) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                        w_L = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)])
                        N_L = 0.0
                        Y_AL_Z_TKE = 0.0
                    else 
                        Y_AL_Z_TKE = y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A_L) .& (y_A_Z_TKE.Z .== Z_L) .& (y_A_Z_TKE.TKE .== TKE)][1]
                        w_L = last(Raw_output_datafile.No_Sequence[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)])
                        T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]
                        α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]
                        E_f = E_f_L(A_0, A_L, A_H, TKE)
                        u_1_L = u_1(E, E_f)
                        u_2_L = u_2(E, E_f)
                        Integrand_L(ε::Float64) = (sum(Neutron_spectrum_k.(ε, α, T)) /w_L)/(4 *sqrt(ε*E_f))
                        N_L = quadgk(Integrand_L, u_1_L, u_2_L)[1]
                    end
                    Denominator += Y_AH_Z_TKE + Y_AL_Z_TKE
                    Numerator += Y_AH_Z_TKE*w_H*N_H + Y_AL_Z_TKE*w_L*N_L
                end
            end
        end
        if Denominator != 0
            push!(n_E.E, E)
            push!(n_E.Value, Numerator/Denominator)
        end
    end
    return n_E
end
function Neutron_spectrum_builder(A_0, Z_0, A_H_min, A_H_max, fragmdomain::Distribution, tkerange, energyrange, Raw_output_datafile::DataFrame, evaporation_cs_type, y_A_Z_TKE::Distribution)
    if evaporation_cs_type == "CONSTANT"
        return Neutron_spectrum_CONSTANT_cs(A_0, Z_0, A_H_min, A_H_max, fragmdomain::Distribution, tkerange, energyrange, Raw_output_datafile::DataFrame, y_A_Z_TKE::Distribution)
    elseif evaporation_cs_type == "VARIABLE"
        return Neutron_spectrum_VARIABLE_cs(A_0, Z_0, A_H_min, A_H_max, fragmdomain::Distribution, tkerange, energyrange, Raw_output_datafile::DataFrame, y_A_Z_TKE::Distribution)
    end
end
#####

n_E = Neutron_spectrum_builder(A₀, Z₀, A_H_min, A_H_max, fragmdomain, tkerange, energyrange, Raw_output_datafile, evaporation_cs_type, y_A_Z_TKE)