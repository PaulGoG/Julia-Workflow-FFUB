#Iteratively compute total neutron spectrum in Laboratory frame from main DSE data
#####
function Neutron_spectrum_k(ε::Float64, T::Float64)
    return ε*exp(-ε/T)/T^2
end
function Neutron_spectrum_k(ε::Float64, α::Float64, T::Float64)
    (ε + α*sqrt(ε))*exp(-ε/T) /((sqrt(T) +α*sqrt(π)/2) *T^(3/2))
end
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
function Define_Integrand(evaporation_cs_type, w, Raw_output_datafile, A, Z, TKE)
    T = Raw_output_datafile.Tₖ[(Raw_output_datafile.A .== A) .& (Raw_output_datafile.Z .== Z) .& (Raw_output_datafile.TKE .== TKE)]
    if evaporation_cs_type == "CONSTANT"
        Neutron_spectrum(ϵ::Float64) = sum(Neutron_spectrum_k.(ϵ, T)) /w
    elseif evaporation_cs_type == "VARIABLE"
        α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A) .& (Raw_output_datafile.Z .== Z) .& (Raw_output_datafile.TKE .== TKE)]
       Neutron_spectrum(ϵ::Float64) = sum(Neutron_spectrum_k.(ϵ, α, T)) /w
    end
end
function Neutron_spectrum_builder(A_0, Z_0, A_H_min, A_H_max, fragmdomain::Distribution, tkerange, energyrange, Raw_output_datafile::DataFrame, evaporation_cs_type, y_A_Z_TKE::Distribution)
    n_E = Neutron_spectrum_struct(Float64[], Float64[])
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
                        if evaporation_cs_type == "CONSTANT"
                            Neutron_spectrum(ϵ::Float64) = sum(Neutron_spectrum_k.(ϵ, T))/w_L
                            Integrand(ϵ::Float64) = Neutron_spectrum(ϵ)/(4 *sqrt(ϵ*E_f))
                            N_H = quadgk(Integrand, u_1_L, u_2_L)[1]
                        elseif evaporation_cs_type == "CONSTANT"
                            α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]
                            Neutron_spectrum(ϵ::Float64) = sum(Neutron_spectrum_k.(ϵ, α, T))/w_L
                            Neutron_spectrum(2.0)
                            Integrand(ϵ::Float64) = Neutron_spectrum(ϵ)/(4 *sqrt(ϵ*E_f))
                            N_H = quadgk(Integrand, u_1_H, u_2_H)[1]
                        end
                        α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A_H) .& (Raw_output_datafile.Z .== Z_H) .& (Raw_output_datafile.TKE .== TKE)]
                        Neutron_spectrum(ϵ::Float64) = sum(Neutron_spectrum_k.(ϵ, α, T))/w_L
                        Neutron_spectrum(2.0)
                        Integrand(ϵ::Float64) = Neutron_spectrum(ϵ)/(4 *sqrt(ϵ*E_f))
                        N_H = quadgk(Integrand, u_1_H, u_2_H)[1]
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
                        if evaporation_cs_type == "CONSTANT"
                            Neutron_spectrum(ϵ::Float64) = sum(Neutron_spectrum_k.(ϵ, T))/w_L
                            Integrand(ϵ::Float64) = Neutron_spectrum(ϵ)/(4 *sqrt(ϵ*E_f))
                            N_L = quadgk(Integrand, u_1_L, u_2_L)[1]
                        elseif evaporation_cs_type == "CONSTANT"
                            α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]
                            Neutron_spectrum(ϵ::Float64) = sum(Neutron_spectrum_k.(ϵ, α, T))/w_L
                            Integrand(ϵ::Float64) = Neutron_spectrum(ϵ)/(4 *sqrt(ϵ*E_f))
                            N_L = quadgk(Integrand, u_1_L, u_2_L)[1]
                        end
                        α = Raw_output_datafile.αₖ[(Raw_output_datafile.A .== A_L) .& (Raw_output_datafile.Z .== Z_L) .& (Raw_output_datafile.TKE .== TKE)]
                        Neutron_spectrum(ϵ::Float64) = sum(Neutron_spectrum_k.(ϵ, α, T))/w_L
                        Integrand(ϵ::Float64) = Neutron_spectrum(ϵ)/(4 *sqrt(ϵ*E_f))
                        N_L = quadgk(Integrand, u_1_L, u_2_L)[1]
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
#####

n_E = Neutron_spectrum_builder(A₀, Z₀, A_H_min, A_H_max, fragmdomain, tkerange, energyrange, Raw_output_datafile, evaporation_cs_type, y_A_Z_TKE)

