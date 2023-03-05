#Function bodies and definitions used in the main part of the DSE
#####
#Load Julia packages
using DataFrames, CSV

#Main struct objects definitions
abstract type AbstractDistribution end
struct Distribution{T1 <: Vector{Int}, T2 <: Vector{Float64}} <: AbstractDistribution
    A::T1
    Z::T1
    TKE::T2
    No_Sequence::T1
    Value
    σ::T2
end
struct Distribution_unidym{T <: Vector{Float64}} <: AbstractDistribution
    Argument
    Value::T
    σ::T
end
struct Neutron_spectrum{T <: Vector{Float64}} <: AbstractDistribution
    ε::T
    Value::T
end

#Define value range for TKE
tkerange = TKE_min:TKE_step:TKE_max

#Input variables corrections according to fission type
if fission_type == "SF"
    #Null neutron incident energy in spontaneous fission
    Eₙ = 0.0
elseif fission_type == "(n,f)"
    #Taking into account compound nucleus formation
    A₀ += 1
end

if evaporation_cs_type == "VARIABLE"
    using Roots
    ħc = 197.3268601
    amu = 931.50176
    aₘ = 1.008665
    const r₀ = 1.2e-1
    const C_α = (π*ħc)^2 /(aₘ*amu)
end

println("*reading data files")
#Read input data files as DataFrames
dmass_excess = CSV.read(mass_excess_filename, DataFrame; delim = mass_excess_delimiter, ignorerepeated = true, header = mass_excess_header, skipto = mass_excess_firstdataline)
println("reading $mass_excess_filename done!")

if density_parameter_type == "GC"
    density_parameter_datafile = CSV.read(density_parameter_filename, DataFrame; delim = density_parameter_delimiter, ignorerepeated = true, header = density_parameter_header, skipto = density_parameter_firstdataline)
    println("reading $density_parameter_filename done!")
elseif density_parameter_type == "BSFG"
    density_parameter_datafile = dmass_excess
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
    isobaric_distribution_datafile = DataFrame(A = NaN)
elseif isobaric_distribution_type == "DATA"
    isobaric_distribution_datafile = CSV.read(isobaric_distribution_filename, DataFrame; delim = isobaric_distribution_delimiter, ignorerepeated = true, header = isobaric_distribution_header, skipto = isobaric_distribution_firstdataline)
    println("reading $isobaric_distribution_filename done!")
end

if txe_partitioning_type == "MSCZ"
    txe_partitioning_datafile = CSV.read(txe_partitioning_filename, DataFrame; delim = txe_partitioning_delimiter, ignorerepeated = true, header = txe_partitioning_header, skipto = txe_partitioning_firstdataline)
    println("reading $txe_partitioning_filename done!")
else
    txe_partitioning_datafile = txe_partitioning_segmentpoints
end

if secondary_outputs == "YES"
    include("secondary_outputs.jl")
end
#Function bodies
#Isobaric charge distribution p(Z,A)
function p_A_Z(Z, Z_p, rms_A)
    return 1/(sqrt(2*π) * rms_A) * exp(-(Z - Z_p)^2 /(2* rms_A^2))
end
#Compute the most probable charge for a given heavy fragment
function Most_probable_charge(A, Z, A_H, ΔZ)
    return A_H*Z/A + ΔZ
end
#Q_value energy in MeV released at fission of (A,Z) nucleus 
function Q_value_released(A_0, Z_0, A_H, Z_H, dm)
    D = dm.D[(dm.A .== A_0) .& (dm.Z .== Z_0)][1]
    σ_D = dm.σ_D[(dm.A .== A_0) .& (dm.Z .== Z_0)][1]
    A_L = A_0 - A_H
    Z_L = Z_0 - Z_H
    if isassigned(dm.D[(dm.A .== A_H) .& (dm.Z .== Z_H)], 1) && isassigned(dm.D[(dm.A .== A_L) .& (dm.Z .== Z_L)], 1)
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
    if isassigned(dm.D[(dm.A .== A_part) .& (dm.Z .== Z_part)], 1) && isassigned(dm.D[(dm.A .== A) .& (dm.Z .== Z)], 1) && isassigned(dm.D[(dm.A .== A - A_part) .& (dm.Z .== Z - Z_part)], 1)
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
#Total Excitation Energy for a given set of data in MeV
function Total_excitation_energy(Q, σ_Q, TKE, σ_TKE, Sₙ, σ_Sₙ, Eₙ)
    TXE = Q - TKE + Sₙ + Eₙ
    σ_TXE = sqrt(σ_Q^2 + σ_Sₙ^2 + σ_TKE^2)
    if TXE > 0
        return TXE, σ_TXE
    else 
        return NaN
    end
end
#Construct vectorized fragmentation domain with p(A,Z) values stored in memory
function Fragmentation_domain(A_0, Z_0, NoZperA, A_H_min, A_H_max, dpAZ)
    fragmdomain = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A_H in A_H_min:A_H_max
        A_L = A_0 - A_H
        if isassigned(dpAZ.A[dpAZ.A .== A_H], 1)
            RMS = dpAZ.rms_A[dpAZ.A .== A_H][1]
            ΔZ = dpAZ.ΔZ_A[dpAZ.A .== A_H][1]
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
#Energy in Fermi Gas regime
function Energy_FermiGas(a::Float64, T::Float64)      
    return a *T^2
end
#Function bodies for neutron spectrum and average neutron energy for a given sequence
if evaporation_cs_type == "CONSTANT"
    function Neutron_spectrum(ε::Float64, T::Float64)
        return ε*exp(-ε/T)/T^2
    end
    function Average_neutron_energy(T::Float64)
        return 2*T
    end
elseif evaporation_cs_type == "VARIABLE"
    function Neutron_spectrum(ε::Float64, α::Float64, T::Float64)
        (ε + α*sqrt(ε))*exp(-ε/T) /((sqrt(T) +α*sqrt(π)/2) *T^(3/2))
    end
    function Average_neutron_energy(α::Float64, T::Float64)
        return T*(2*sqrt(T) + α*3*sqrt(π)/4) /(sqrt(T) + α*sqrt(π)/2)
    end
end
#Processing raw DSE eq output
function Process_main_output(DSE_eq_output, evaporation_cs_type)
    Tₖ, aₖ = DSE_eq_output[1], DSE_eq_output[2]
    if evaporation_cs_type .== "CONSTANT"
        Datafile = DataFrame(
            A = Tₖ.A, 
            Z = Tₖ.Z, 
            TKE = Tₖ.TKE, 
            No_Sequence = Tₖ.No_Sequence, 
            Tₖ = Tₖ.Value, 
            aₖ = aₖ
        )
    elseif evaporation_cs_type .== "VARIABLE"
        αₖ = DSE_eq_output[3]
        Datafile = DataFrame(
            A = Tₖ.A, 
            Z = Tₖ.Z, 
            TKE = Tₖ.TKE, 
            No_Sequence = Tₖ.No_Sequence, 
            Tₖ = Tₖ.Value, 
            aₖ = aₖ,
            αₖ = αₖ
        )
    end
    return Datafile
end
#Writing main output file
function Write_seq_output(A_0, Z_0, A_H_min, A_H_max, No_ZperA, Eₙ, tkerange, fragmdomain, E_excitation, Processed_raw_output, density_parameter_type, density_parameter_datafile, evaporation_cs_type, mass_excess_filename, txe_partitioning_type, dm)
    if !isdir("output_data/")
        mkdir("output_data/")
    end
    Sₙ = Separation_energy(1, 0, A_0, Z_0, dm)
    open("output_data/main_DSE.OUT", "w") do file
        write(file, "DSE main output file for the following input data:\n")
        write(file,"(A₀ = $A_0, Z₀ = $Z_0), $fission_type fission, $No_ZperA Z per A, mass excess file - $mass_excess_filename\n")
        write(file,"Heavy Fragment mass number range from $A_H_min to $A_H_max\n")
        write(file,"TKE range ∈ $tkerange\n")
        write(file,"TXE partitioning method - $txe_partitioning_type\n")
        write(file, "\n===========================================================================================================================================================\n")
        for A in unique(Processed_raw_output.A)
            for Z in unique(Processed_raw_output.Z[Processed_raw_output.A .== A])
                P_Z_A = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
                Q = Q_value_released(A_0, Z_0, A, Z, dm)
                a = density_parameter(density_parameter_type, A, Z, density_parameter_datafile)
                S = Separation_energy(1, 0, A, Z, dm)[1]
                write(file, "Fission fragment: A = $A Z = $Z p(Z,A) = $P_Z_A Q = $(Q[1]) a = $a Sₙ = $S\n\n")
                for TKE in unique(Processed_raw_output.TKE[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z)])
                    TXE = Total_excitation_energy(Q[1], Q[2], TKE, 0.0, Sₙ[1], Sₙ[2], Eₙ)[1]
                    E_excit = E_excitation.Value[(E_excitation.A .== A) .& (E_excitation.Z .== Z) .& (E_excitation.TKE .== TKE)][1]
                    write(file, "TKE = $TKE TXE = $TXE E* = $E_excit\n")
                    for k in unique(Processed_raw_output.No_Sequence[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z) .& (Processed_raw_output.TKE .== TKE)])
                        write(file, "   *Emission sequence $k:\n   ")
                        T_k = Processed_raw_output.Tₖ[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                        write(file, "T = $T_k ")
                        a_k = Processed_raw_output.aₖ[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                        write(file, "a = $a_k ")
                        Eᵣ_k = Energy_FermiGas(a_k, T_k)
                        write(file, "Eʳ = $Eᵣ_k ")
                        S_k = Separation_energy(1, 0, A-k, Z, dm)[1]
                        write(file, "Sₙ = $S_k ")
                        if evaporation_cs_type == "CONSTANT"
                            avgε_k = Average_neutron_energy(T_k)
                            write(file, "<ε> = $avgε_k\n")
                        elseif evaporation_cs_type == "VARIABLE"
                            α_k = Processed_raw_output.αₖ[(Processed_raw_output.A .== A) .& (Processed_raw_output.Z .== Z) .& (Processed_raw_output.TKE .== TKE) .& (Processed_raw_output.No_Sequence .== k)][1]
                            avgε_k = Average_neutron_energy(α_k, T_k)
                            write(file, "<ε> = $avgε_k\n")
                        end
                    end
                    write(file, '\n')
                end
                write(file, "-----------------------------------------------------------------------------------------------------------------------------------------------------------\n")
            end
        end
    end
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
                push!(ν.Value, last(output_df_A_Z_TKE_NoSequence.No_Sequence[(output_df_A_Z_TKE_NoSequence.A .== A) .& (output_df_A_Z_TKE_NoSequence.Z .== Z) .& (output_df_A_Z_TKE_NoSequence.TKE .== TKE)]))
            end
        end
    end
    return ν
end
#Average raw output data over emission sequences
function SeqAvg_A_Z_TKE(output_df_A_Z_TKE_NoSequence_Value::DataFrame)
    avg = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A in unique(output_df_A_Z_TKE_NoSequence_Value.A)
        for Z in unique(output_df_A_Z_TKE_NoSequence_Value.Z[output_df_A_Z_TKE_NoSequence_Value.A .== A])
            for TKE in unique(output_df_A_Z_TKE_NoSequence_Value.TKE[(output_df_A_Z_TKE_NoSequence_Value.A .== A) .& (output_df_A_Z_TKE_NoSequence_Value.Z .== Z)])
                push!(avg.A, A)
                push!(avg.Z, Z)
                push!(avg.TKE, TKE)
                n = last(output_df_A_Z_TKE_NoSequence_Value.No_Sequence[(output_df_A_Z_TKE_NoSequence_Value.A .== A) .& (output_df_A_Z_TKE_NoSequence_Value.Z .== Z) .& (output_df_A_Z_TKE_NoSequence_Value.TKE .== TKE)])
                push!(avg.Value, sum(output_df_A_Z_TKE_NoSequence_Value.Value[(output_df_A_Z_TKE_NoSequence_Value.A .== A) .& (output_df_A_Z_TKE_NoSequence_Value.Z .== Z) .& (output_df_A_Z_TKE_NoSequence_Value.TKE .== TKE)])/n)
            end
        end
    end
    return avg
end