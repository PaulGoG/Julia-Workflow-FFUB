#=
Function bodies for solving the DSE conservation equations coresponding to 
constant and variable neutron evaporation cross section types.
=#
#####
function DSE_equation_solver_CONSTANT_cs(fragmdomain, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    Tₖ = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    aₖ = Float64[]
    for A in unique(fragmdomain.A)
        for Z in fragmdomain.Z[fragmdomain.A .== A]
            Sₙ = Separation_energy(1, 0, A, Z, dm)[1]
            a_1 = density_parameter(density_parameter_type, A - 1, Z, density_parameter_datafile)
            if !isnan(a_1)
                for TKE in E_excitation.TKE[(E_excitation.A .== A) .& (E_excitation.Z .== Z)]
                    Eᵣ_k_last = E_excitation.Value[(E_excitation.A .== A) .& (E_excitation.Z .== Z) .& (E_excitation.TKE .== TKE)][1]
                    Sₙ_k_last = Sₙ
                    a_k = a_1
                    k = 1
                    while Eᵣ_k_last > Sₙ_k_last
                        T_k = (sqrt(1 + a_k * (Eᵣ_k_last - Sₙ_k_last)) - 1)/a_k
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, T_k)
                        push!(aₖ, a_k)
                        push!(Tₖ.NoSeq, k)
                        #Advance the sequence one step forward to be verified by the while loop
                        Eᵣ_k_last = a_k *T_k^2
                        Sₙ_k_last = Separation_energy(1, 0, A - k, Z, dm)[1]
                        k += 1
                        a_k = density_parameter(density_parameter_type, A - k, Z, density_parameter_datafile)
                    end
                end  
            end          
        end
    end
    return Tₖ, aₖ
end
#VARIABLE σₙ functions
#Parametrization for the force function S₀ of the s-wave neutron
function Force_function_S₀(A)
    if A <= 70
        return 7e-5
    elseif A > 70 && A <= 86
        return 1e-4 * (A*1.875e-2 - 6.125e-1)
    elseif A > 86 && A <= 111
        return 1e-4 
    elseif A > 111 && A <= 121
        return 1e-4 * (-A*2.857e-2 + 4.1714)
    elseif A > 121 && A <= 140
        return 1e-4 
    elseif A > 140 && A <= 144
        return 1e-4 * (A*7.5e-2 - 9.8)
    elseif A > 144
        return 1e-4
    end
end
#Solves the transcendental equation for given sequence
function Solve_transcendental_eq(Eᵣ_k_last, Sₙ_k_last, a_k, A, k)
    σ₀ = π*r₀^2 *(A - k)^(2/3)
    S₀ = Force_function_S₀(A - k)
    αₖ = 10*C_α*S₀/σ₀
    f(Tₖ) = a_k*Tₖ^2 + Tₖ*(2*sqrt(Tₖ) + αₖ*3*sqrt(π)/4)/(sqrt(Tₖ) + αₖ*sqrt(π)/2) + (Sₙ_k_last - Eᵣ_k_last)
    T_k = find_zero(f, 1.0)
    return T_k, αₖ
end
function DSE_equation_solver_VARIABLE_cs(fragmdomain, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    Tₖ = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    aₖ = Float64[]
    αₖ = Float64[]
    for A in first(fragmdomain.A):last(fragmdomain.A)
        for Z in fragmdomain.Z[fragmdomain.A .== A]
            Sₙ = Separation_energy(1, 0, A, Z, dm)[1]
            a_1 = density_parameter(density_parameter_type, A - 1, Z, density_parameter_datafile)
            if !isnan(a_1)
                for TKE in E_excitation.TKE[(E_excitation.A .== A) .& (E_excitation.Z .== Z)]
                    Eᵣ_k_last = E_excitation.Value[(E_excitation.A .== A) .& (E_excitation.Z .== Z) .& (E_excitation.TKE .== TKE)][1]
                    Sₙ_k_last = Sₙ
                    a_k = a_1
                    k = 1
                    while Eᵣ_k_last > Sₙ_k_last
                        T_k, α_k = Solve_transcendental_eq(Eᵣ_k_last, Sₙ_k_last, a_k, A, k)
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, T_k)
                        push!(aₖ, a_k)
                        push!(αₖ, α_k)
                        push!(Tₖ.NoSeq, k)
                        #Advance the sequence one step forward to be verified by the while loop
                        Eᵣ_k_last = a_k *T_k^2
                        Sₙ_k_last = Separation_energy(1, 0, A - k, Z, dm)[1]
                        k += 1
                        a_k = density_parameter(density_parameter_type, A - k, Z, density_parameter_datafile)
                    end
                end  
            end          
        end
    end
    return Tₖ, aₖ, αₖ
end
function DSE_equation_solver(evaporation_cs_type, fragmdomain, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    if evaporation_cs_type == "CONSTANT"
        DSE_Output = DSE_equation_solver_CONSTANT_cs(fragmdomain, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    elseif evaporation_cs_type == "VARIABLE"
        DSE_Output = DSE_equation_solver_VARIABLE_cs(fragmdomain, E_excitation, tkerange, density_parameter_type, density_parameter_datafile, dm)
    end
    return DSE_Output
end
