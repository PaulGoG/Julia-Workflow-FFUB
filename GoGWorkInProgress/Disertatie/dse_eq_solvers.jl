#=
Function bodies for solving the DSE energy conservation equations coresponding to 
constant and variable neutron evaporation cross section types.
=#
#####
#CONSTANT σₙ function
function DSE_equation_solver_CONSTANT_cs(E_excitation::Distribution, density_parameter_type, density_parameter_data, dm::DataFrame)
    Tₖ = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    aₖ = Float64[]
    for A in unique(E_excitation.A)
        for Z in unique(E_excitation.Z[E_excitation.A .== A])
            Sₙ = Separation_energy(1, 0, A, Z, dm)[1]
            a_1 = density_parameter(density_parameter_type, A - 1, Z, density_parameter_data)
            if !isnan(a_1)
                for TKE in E_excitation.TKE[(E_excitation.A .== A) .& (E_excitation.Z .== Z)]
                    Eᵣ_k_last = E_excitation.Value[(E_excitation.A .== A) .& (E_excitation.Z .== Z) .& (E_excitation.TKE .== TKE)][1]
                    Sₙ_k_last = Sₙ
                    a_k = a_1
                    k = 1
                    while Eᵣ_k_last > Sₙ_k_last
                        T_k = (sqrt(1 + a_k *(Eᵣ_k_last - Sₙ_k_last)) - 1) /a_k
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, T_k)
                        push!(aₖ, a_k)
                        push!(Tₖ.No_Sequence, k)
                        #Advance the sequence one step forward to be verified by the while loop
                        Eᵣ_k_last = Energy_FermiGas(a_k, T_k)
                        Sₙ_k_last = Separation_energy(1, 0, A - k, Z, dm)[1]
                        k += 1
                        a_k = density_parameter(density_parameter_type, A - k, Z, density_parameter_data)
                    end
                    if k == 1
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, NaN)
                        push!(aₖ, a_k)
                        push!(Tₖ.No_Sequence, 0)
                    end
                end  
            end          
        end
    end
    return Tₖ, aₖ
end
#VARIABLE σₙ functions
#Parametrization for the force function S₀ of the s-wave neutron
function Neutron_strength_function(A)
    if A <= 23
        return 2.5e-5
    elseif A > 23 && A <= 55
        return 1.296875e-5 *A - 0.00027328125
    elseif A > 55 && A <= 69
        return -2.142857142857143e-5 *A + 0.0016185714285714289
    elseif A > 69 && A <= 79
        return 0.00014
    elseif A > 79 && A <= 90
        return -7.818181818181818e-6 *A + 0.0007576363636363635
    elseif A > 90 && A <= 110
        return 5.4e-5
    elseif A > 110 && A <= 120
        return -4.4e-6 *A + 0.000538 
    elseif A > 120 && A <= 135
        return 6.666666666666667e-6 *A - 0.00079
    elseif A > 135 && A <= 146
        return 2.4545454545454552e-5 *A - 0.0032036363636363647
    elseif A > 146 && A <= 152
        return 0.00038
    elseif A > 152 && A <= 158
        return -3.833333333333334e-5 *A + 0.006206666666666668
    elseif A > 158 && A <= 173
        return 0.00015
    elseif A > 173 && A <= 186
        return 1.0769230769230771e-5 *A - 0.0017130769230769235
    elseif A > 186 && A <= 208
        return -8.636363636363637e-6 *A + 0.0018963636363636366
    elseif A > 208
        return 1e-4
    end
end
#Solves the transcendental equation for given sequence
function Solve_transcendental_eq(Eᵣ_k_last, Sₙ_k_last, a_k, A, k)
    σ₀ = π *(r₀ *0.1)^2 *(A - k)^(2/3)
    S₀ = Neutron_strength_function(A - k)
    αₖ = 10 *C_α *S₀ /σ₀
    f(Tₖ) = a_k*Tₖ^2 + Tₖ*(2*sqrt(Tₖ) + αₖ*3*sqrt(π)/4)/(sqrt(Tₖ) + αₖ*sqrt(π)/2) + (Sₙ_k_last - Eᵣ_k_last)
    T_k = find_zero(f, 1.0)
    return T_k, αₖ
end
function DSE_equation_solver_VARIABLE_cs(E_excitation::Distribution, density_parameter_type, density_parameter_data, dm::DataFrame)
    Tₖ = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    aₖ = Float64[]
    αₖ = Float64[]
    for A in unique(E_excitation.A)
        for Z in unique(E_excitation.Z[E_excitation.A .== A])
            Sₙ = Separation_energy(1, 0, A, Z, dm)[1]
            a_1 = density_parameter(density_parameter_type, A - 1, Z, density_parameter_data)
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
                        push!(Tₖ.No_Sequence, k)
                        #Advance the sequence one step forward to be verified by the while loop
                        Eᵣ_k_last = Energy_FermiGas(a_k, T_k)
                        Sₙ_k_last = Separation_energy(1, 0, A - k, Z, dm)[1]
                        k += 1
                        a_k = density_parameter(density_parameter_type, A - k, Z, density_parameter_data)
                    end
                    if k == 1
                        push!(Tₖ.A, A)
                        push!(Tₖ.Z, Z)
                        push!(Tₖ.TKE, TKE)
                        push!(Tₖ.Value, NaN)
                        push!(aₖ, a_k)
                        push!(αₖ, NaN)
                        push!(Tₖ.No_Sequence, 0)
                    end
                end  
            end          
        end
    end
    return Tₖ, aₖ, αₖ
end
function DSE_equation_solver(evaporation_cs_type, E_excitation, density_parameter_type, density_parameter_data, dm)
    println("*solving DSE energy conservation equations")
    if evaporation_cs_type == "CONSTANT"
        return DSE_equation_solver_CONSTANT_cs(E_excitation, density_parameter_type, density_parameter_data, dm)
    elseif evaporation_cs_type == "VARIABLE"
        return DSE_equation_solver_VARIABLE_cs(E_excitation, density_parameter_type, density_parameter_data, dm)
    end
end
