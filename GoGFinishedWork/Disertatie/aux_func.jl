#Define value range for TKE
tkerange = TKE_min:TKE_step:TKE_max

#Input variables corrections according to fission type
if fission_type == "SF"
    #Null neutron incident energy in spontaneous fission
    E_incident = 0.0
elseif fission_type == "(n,f)"
    #Take into account compound nucleus formation
    A₀ += 1
end

#Define value ranges for fragmentation regions
A_H_range = A_H_min:A_H_max
A_L_range = (A₀ - A_H_max):(A₀ - A_H_min)
A_range = (A₀ - A_H_max):A_H_max

#Read data and initialise constants where needed

if evaporation_cs_type == "VARIABLE"
    using Roots
    ħc = 197.3268601
    amu = 931.50176
    aₘ = 1.008665
    const r₀ = 1.2
    const C_α = (π*ħc)^2 /(aₘ*amu)
end

if neutron_spectrum
    using QuadGK, Trapz
    struct Neutron_spectrum{T <: Vector{Float64}} <: AbstractDistribution
        E::T
        Value::T
    end
    energyrange = E_min:E_step:E_max
end

dmass_excess = CSV.read(mass_excess_filename, DataFrame; 
delim = mass_excess_delimiter, ignorerepeated = true, header = mass_excess_header, skipto = mass_excess_firstdataline)
println("*reading $mass_excess_filename done!")

if density_parameter_type == "GC"
    density_parameter_data = CSV.read(density_parameter_filename, DataFrame; 
    delim = density_parameter_delimiter, ignorerepeated = true, header = density_parameter_header, skipto = density_parameter_firstdataline)
    println("*reading $density_parameter_filename done!")
elseif density_parameter_type == "BSFG"
    density_parameter_data = dmass_excess
    const Dᵖ = dmass_excess.D[(dmass_excess.A .== 1) .& (dmass_excess.Z .== 1)][1]*1e-3
    const Dⁿ = dmass_excess.D[(dmass_excess.A .== 1) .& (dmass_excess.Z .== 0)][1]*1e-3
    const a_v = 15.65
    const a_s = 17.63
    const a_c = 0.864/1.233
    const a_sim_1 = 27.72
    const a_sim_2 = 25.6
    const p_eb_1 = 1.99e-1
    const p_eb_2 = 9.6e-3
    const p_eb_3 = 8.69e-1
end

if isobaric_distribution_type == "MEAN_VALUES"
    isobaric_distribution_data = DataFrame(A = NaN)
elseif isobaric_distribution_type == "DATA"
    isobaric_distribution_data = CSV.read(isobaric_distribution_filename, DataFrame; 
    delim = isobaric_distribution_delimiter, ignorerepeated = true, header = isobaric_distribution_header, skipto = isobaric_distribution_firstdataline)
    println("*reading $isobaric_distribution_filename done!")
end

if txe_partitioning_type == "MSCZ"
    txe_partitioning_data = CSV.read(txe_partitioning_filename, DataFrame; 
    delim = txe_partitioning_delimiter, ignorerepeated = true, header = txe_partitioning_header, skipto = txe_partitioning_firstdataline)
    println("*reading $txe_partitioning_filename done!")
else
    txe_partitioning_data = txe_partitioning_segmentpoints
end

if secondary_outputs
    Yield_data = CSV.read(yield_distribution_filename, DataFrame; 
    delim = yield_distribution_delimiter, ignorerepeated = true, header = yield_distribution_header, skipto = yield_distribution_firstdataline)
    println("*reading $yield_distribution_filename done!")
end

#Revert relative PATH to project root folder
cd(@__DIR__)

if !isdir("$(file_output_identifier)_output_data/")
    mkdir("$(file_output_identifier)_output_data/")
end

if generate_plots
    using Plots, LaTeXStrings
    plots_resolution = aspect_ratio .* resolution_scale
    if !isdir("$(file_output_identifier)_plots/")
        mkdir("$(file_output_identifier)_plots/")
    end
end

#Define main functions

#Most probable charge for a given heavy fragment
function Most_probable_charge(A_0, Z_0, A_H, ΔZ)
    return ΔZ + A_H *Z_0 /A_0
end
#Isobaric charge distribution p(Z,A)
function p_A_Z(Z, Z_p, rms_A)
    return exp(-(Z - Z_p)^2 /(2 *rms_A^2)) /(sqrt(2*π) *rms_A)
end
#Q_value energy in MeV released at fission of (A,Z) nucleus 
function Q_value_released(A_0, Z_0, A_H, Z_H, dm)
    D = dm.D[(dm.A .== A_0) .& (dm.Z .== Z_0)][1]
    σ_D = dm.σ_D[(dm.A .== A_0) .& (dm.Z .== Z_0)][1]
    A_L = A_0 - A_H
    Z_L = Z_0 - Z_H
    if isassigned(dm.D[(dm.A .== A_H) .& (dm.Z .== Z_H)], 1) && 
        isassigned(dm.D[(dm.A .== A_L) .& (dm.Z .== Z_L)], 1)

        D_H = dm.D[(dm.A .== A_H) .& (dm.Z .== Z_H)][1]
        σ_D_H = dm.σ_D[(dm.A .== A_H) .& (dm.Z .== Z_H)][1]
        D_L = dm.D[(dm.A .== A_L) .& (dm.Z .== Z_L)][1]
        σ_D_L = dm.σ_D[(dm.A .== A_L) .& (dm.Z .== Z_L)][1]
        Q = (D - (D_H + D_L)) *1e-3
        σ = sqrt(σ_D^2 + σ_D_H^2 + σ_D_L^2) *1e-3
        if Q > 0
            return Q, σ
        else
            return NaN
        end
    else
        return NaN
    end
end
#Separation energy of particle (A_part,Z_part) from nucleus (A,Z) in MeV
function Separation_energy(A_part, Z_part, A, Z, dm)
    if isassigned(dm.D[(dm.A .== A_part) .& (dm.Z .== Z_part)], 1) && 
        isassigned(dm.D[(dm.A .== A) .& (dm.Z .== Z)], 1) && 
        isassigned(dm.D[(dm.A .== A - A_part) .& (dm.Z .== Z - Z_part)], 1)

        D_part = dm.D[(dm.A .== A_part) .& (dm.Z .== Z_part)][1]
        σ_part = dm.σ_D[(dm.A .== A_part) .& (dm.Z .== Z_part)][1]
        D = dm.D[(dm.A .== A - A_part) .& (dm.Z .== Z - Z_part)][1]
        σ_D = dm.σ_D[(dm.A .== A - A_part) .& (dm.Z .== Z - Z_part)][1]    
        S = (D + D_part - dm.D[(dm.A .== A) .& (dm.Z .== Z)][1]) *1e-3
        σ = sqrt(σ_part^2 + σ_D^2 + (dm.σ_D[(dm.A .== A) .& (dm.Z .== Z)][1])^2) *1e-3
        return S, σ
    else 
        return NaN
    end 
end
#Excitation energy of compound nucleus
function Compound_nucleus_energy(fission_type, A_0, Z_0, E_incident, dm)
    if fission_type == "SF"
        return 0.0, 0.0
    elseif fission_type == "(n,f)" 
        S = Separation_energy(1, 0, A_0, Z_0, dm)
        return S[1] + E_incident, S[2]
    end
end
#Total Excitation Energy for a given set of data in MeV
function Total_excitation_energy(Q, σ_Q, TKE, σ_TKE, ε_CN, σ_ε_CN)
    TXE = Q - TKE + ε_CN
    σ_TXE = sqrt(σ_Q^2 + σ_ε_CN^2 + σ_TKE^2)
    if TXE > 0
        return TXE, σ_TXE
    else 
        return NaN
    end
end
#Construct vectorized fragmentation domain with p(A,Z) values stored in memory
function Fragmentation_domain(A_0, Z_0, NoZperA, A_H_range, pAZ_data)
    println("*building fragmentation domain")
    fragmdomain = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A_H in A_H_range
        A_L = A_0 - A_H
        if isassigned(pAZ_data.A[pAZ_data.A .== A_H], 1)
            RMS = pAZ_data.rms_A[pAZ_data.A .== A_H][1]
            ΔZ = pAZ_data.ΔZ_A[pAZ_data.A .== A_H][1]
        else
            RMS = 0.6
            ΔZ = -0.5
        end
        Zₚ = Most_probable_charge(A_0, Z_0, A_H, ΔZ)
        Z_H_min = Int(round(Zₚ) - (NoZperA - 1)/2)
        Z_H_max = Z_H_min + NoZperA - 1
        for Z_H in Z_H_min:Z_H_max
            push!(fragmdomain.A, A_H)
            push!(fragmdomain.Z, Z_H)
            push!(fragmdomain.Value, p_A_Z(Z_H, Zₚ, RMS))
            if A_L != A_H_min
                Z_L = Z_0 - Z_H
                push!(fragmdomain.A, A_L)
                push!(fragmdomain.Z, Z_L)
                push!(fragmdomain.Value, p_A_Z(Z_L, Z_0 - Zₚ, RMS))
            end
        end
    end
    #Sort by mass number in ascending order
    aux_A = sort(fragmdomain.A)
    aux_Z = [sort(fragmdomain.Z[fragmdomain.A .== A]) for A in unique(aux_A)]
    aux_Z = reduce(vcat, aux_Z)
    aux_Value = zeros(length(fragmdomain.Value))
    for index in eachindex(aux_A)
        aux_Value[index] = fragmdomain.Value[(fragmdomain.A .== aux_A[index]) .& (fragmdomain.Z .== aux_Z[index])][1]
    end
    for index in eachindex(aux_A)
        fragmdomain.A[index] = aux_A[index]
        fragmdomain.Z[index] = aux_Z[index]
        fragmdomain.Value[index] = aux_Value[index]
    end
    return fragmdomain
end
#Energy in Fermi Gas model
function Energy_FermiGas(a::Float64, T::Float64)      
    return a *T^2
end
#Average neutron energy for a given sequence
function Average_neutron_energy(T::Float64)
    return 2*T
end
function Average_neutron_energy(α::Float64, T::Float64)
    return T *(2*sqrt(T) + α*3*sqrt(π)/4) /(sqrt(T) + α*sqrt(π)/2)
end
#Processing output of the main DSE equations
function Process_main_output(DSE_eq_output, evaporation_cs_type)
    println("*processing DSE equations primary output")
    Tₖ, aₖ = DSE_eq_output[1], DSE_eq_output[2]
    if evaporation_cs_type == "CONSTANT"
        Data = DataFrame(
            A = Tₖ.A, 
            Z = Tₖ.Z, 
            TKE = Tₖ.TKE, 
            No_Sequence = Tₖ.No_Sequence, 
            Tₖ = Tₖ.Value, 
            aₖ = aₖ,
            Avg_εₖ = Average_neutron_energy.(Tₖ.Value)
        )
    elseif evaporation_cs_type == "VARIABLE"
        αₖ = DSE_eq_output[3]
        Data = DataFrame(
            A = Tₖ.A, 
            Z = Tₖ.Z, 
            TKE = Tₖ.TKE, 
            No_Sequence = Tₖ.No_Sequence, 
            Tₖ = Tₖ.Value, 
            aₖ = aₖ,
            αₖ = αₖ,
            Avg_εₖ = Average_neutron_energy.(αₖ, Tₖ.Value)
        )
    end
    return Data
end
#Sorts any q(A,Z,TKE) in ascending order of A, then Z, then TKE
function Sort_q_A_Z_TKE(q_A_Z_TKE, fragmdomain)
    aux_A = copy(q_A_Z_TKE.A)
    aux_Z = copy(q_A_Z_TKE.Z)
    aux_TKE = zeros(length(q_A_Z_TKE.TKE))
    aux_Value = zeros(length(q_A_Z_TKE.Value))
    aux_index = 1
    for index_fragmdomain in eachindex(fragmdomain.A)
        for TKE in sort(unique(q_A_Z_TKE.TKE[(q_A_Z_TKE.A .== fragmdomain.A[index_fragmdomain]) .& (q_A_Z_TKE.Z .== fragmdomain.Z[index_fragmdomain])]))
            aux_A[aux_index] = fragmdomain.A[index_fragmdomain]
            aux_Z[aux_index] = fragmdomain.Z[index_fragmdomain]
            aux_TKE[aux_index] = TKE
            aux_Value[aux_index] = q_A_Z_TKE.Value[(q_A_Z_TKE.A .== fragmdomain.A[index_fragmdomain]) .& (q_A_Z_TKE.Z .== fragmdomain.Z[index_fragmdomain]) .& (q_A_Z_TKE.TKE .== TKE)][1]
            aux_index += 1
        end
    end
    for index in eachindex(aux_Value)
        q_A_Z_TKE.A[index] = aux_A[index]
        q_A_Z_TKE.Z[index] = aux_Z[index]
        q_A_Z_TKE.TKE[index] = aux_TKE[index]
        q_A_Z_TKE.Value[index] = aux_Value[index]
    end
end
#Maximum number of sequences from raw output data
function Maximum_sequences_A_Z_TKE(output_df_A_Z_TKE_NoSequence::DataFrame)
    sequences = Distribution(Int[], Int[], Float64[], Int[], Int[], Float64[])
    for A in unique(output_df_A_Z_TKE_NoSequence.A)
        for Z in unique(output_df_A_Z_TKE_NoSequence.Z[output_df_A_Z_TKE_NoSequence.A .== A])
            for TKE in unique(output_df_A_Z_TKE_NoSequence.TKE[(output_df_A_Z_TKE_NoSequence.A .== A) .& (output_df_A_Z_TKE_NoSequence.Z .== Z)])
                val = last(output_df_A_Z_TKE_NoSequence.No_Sequence[(output_df_A_Z_TKE_NoSequence.A .== A) .& (output_df_A_Z_TKE_NoSequence.Z .== Z) .& (output_df_A_Z_TKE_NoSequence.TKE .== TKE)])
                push!(sequences.A, A)
                push!(sequences.Z, Z)
                push!(sequences.TKE, TKE)
                push!(sequences.Value, val)
            end
        end
    end
    return sequences
end
#Neutron multiplicity from raw output data
function Neutron_multiplicity_A_Z_TKE(output_df_A_Z_TKE_NoSequence::DataFrame)
    ν = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(output_df_A_Z_TKE_NoSequence.A)
        for Z in unique(output_df_A_Z_TKE_NoSequence.Z[output_df_A_Z_TKE_NoSequence.A .== A])
            for TKE in unique(output_df_A_Z_TKE_NoSequence.TKE[(output_df_A_Z_TKE_NoSequence.A .== A) .& (output_df_A_Z_TKE_NoSequence.Z .== Z)])
                push!(ν.A, A)
                push!(ν.Z, Z)
                push!(ν.TKE, TKE)
                val = last(output_df_A_Z_TKE_NoSequence.No_Sequence[(output_df_A_Z_TKE_NoSequence.A .== A) .& (output_df_A_Z_TKE_NoSequence.Z .== Z) .& (output_df_A_Z_TKE_NoSequence.TKE .== TKE)])
                if val == 0
                    push!(ν.Value, val)
                else
                    push!(ν.Value, (val + 1) /2)
                end
            end
        end
    end
    return ν
end
function Neutron_multiplicity_Pair_A_Z_TKE(A_0, Z_0, A_H_range, ν_A_Z_TKE::Distribution)
    ν_Pair = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in A_H_range
        for Z in unique(ν_A_Z_TKE.Z[ν_A_Z_TKE.A .== A])
            for TKE in unique(ν_A_Z_TKE.TKE[(ν_A_Z_TKE.A .== A) .& (ν_A_Z_TKE.Z .== Z)])
                if !isassigned(ν_Pair.Value[(ν_Pair.A .== A) .& (ν_Pair.Z .== Z) .& (ν_Pair.TKE .== TKE)], 1)
                    nu_L = ν_A_Z_TKE.Value[(ν_A_Z_TKE.A .== A_0 - A) .& (ν_A_Z_TKE.Z .== Z_0 - Z) .& (ν_A_Z_TKE.TKE .== TKE)][1]
                    nu_H = ν_A_Z_TKE.Value[(ν_A_Z_TKE.A .== A) .& (ν_A_Z_TKE.Z .== Z) .& (ν_A_Z_TKE.TKE .== TKE)][1]
                    nu_Pair = nu_L + nu_H
                    push!(ν_Pair.A, A)
                    push!(ν_Pair.Z, Z)
                    push!(ν_Pair.TKE, TKE)
                    push!(ν_Pair.Value, nu_Pair)
                end
            end
        end
    end
    return ν_Pair
end
#Average raw output data over valid emission sequences
function SeqAvg_A_Z_TKE(output_df_A_Z_TKE_NoSequence_Value::DataFrame)
    avg = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(output_df_A_Z_TKE_NoSequence_Value.A)
        for Z in unique(output_df_A_Z_TKE_NoSequence_Value.Z[output_df_A_Z_TKE_NoSequence_Value.A .== A])
            for TKE in unique(output_df_A_Z_TKE_NoSequence_Value.TKE[(output_df_A_Z_TKE_NoSequence_Value.A .== A) .& (output_df_A_Z_TKE_NoSequence_Value.Z .== Z)])
                if isassigned(output_df_A_Z_TKE_NoSequence_Value.No_Sequence[(output_df_A_Z_TKE_NoSequence_Value.A .== A) .& (output_df_A_Z_TKE_NoSequence_Value.Z .== Z) .& (output_df_A_Z_TKE_NoSequence_Value.TKE .== TKE) .& (output_df_A_Z_TKE_NoSequence_Value.Value .>= 0)], 1)
                    n = last(output_df_A_Z_TKE_NoSequence_Value.No_Sequence[(output_df_A_Z_TKE_NoSequence_Value.A .== A) .& (output_df_A_Z_TKE_NoSequence_Value.Z .== Z) .& (output_df_A_Z_TKE_NoSequence_Value.TKE .== TKE) .& (output_df_A_Z_TKE_NoSequence_Value.Value .>= 0)])
                    val = sum(filter(!isnan, output_df_A_Z_TKE_NoSequence_Value.Value[(output_df_A_Z_TKE_NoSequence_Value.A .== A) .& (output_df_A_Z_TKE_NoSequence_Value.Z .== Z) .& (output_df_A_Z_TKE_NoSequence_Value.TKE .== TKE)]))
                    push!(avg.A, A)
                    push!(avg.Z, Z)
                    push!(avg.TKE, TKE)
                    push!(avg.Value, val /n)
                end
            end
        end
    end
    return avg
end