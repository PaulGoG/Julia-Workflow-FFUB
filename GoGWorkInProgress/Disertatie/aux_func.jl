#Function bodies and definitions used in the main program
#####
#Load Julia packages
using DataFrames, CSV

#Main struct objects definitions
abstract type AbstractDistribution end
struct Distribution{T1 <: Vector{Int}, T2 <: Vector{Float64}} <: AbstractDistribution
    A::T1
    Z::T1
    TKE::T2
    NoSeq::T1
    Value::T2
    σ::T2
end
struct Distribution_unidym{T <: Vector{Float64}} <: AbstractDistribution
    Index::T
    Value::T
    σ::T
end
struct Energy_spectrum{T <: Vector{Float64}} <: AbstractDistribution
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

if txe_partitioning_type == "RT"
    println("Provide input value for R_T = T_L/T_H for fission fragments at total acceleration in Fermi Gas regime...")
    println(' ')
    print("R_T = ")
    const R_T² = parse(Float64, readline())^2
    println(' ')
    txe_partitioning_datafile = DataFrame
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
elseif txe_partitioning_type == "PARAM"
    txe_partitioning_datafile = CSV.read(txe_partitioning_filename, DataFrame; delim = txe_partitioning_delimiter, ignorerepeated = true, header = txe_partitioning_header, skipto = txe_partitioning_firstdataline)
end
println("reading $txe_partitioning_filename done!")
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
function Q_value_released(A, Z, A_H, Z_H, dm)
    D = dm.D[(dm.A .== A) .& (dm.Z .== Z)][1]
    σ_D = dm.σ_D[(dm.A .== A) .& (dm.Z .== Z)][1]
    A_L = A - A_H
    Z_L = Z - Z_H
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
function Fragmentation_domain(A, Z, NoZperA, A_H_min, A_H_max, dpAZ)
    fragmdomain = Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    for A_H in A_H_min:A_H_max
        A_L = A - A_H
        if isassigned(dpAZ.A[dpAZ.A .== A_H], 1)
            RMS = dpAZ.rms_A[dpAZ.A .== A_H][1]
            ΔZ = dpAZ.ΔZ_A[dpAZ.A .== A_H][1]
        else
            RMS = 0.6
            ΔZ = -0.5
        end
        Zₚ = Most_probable_charge(A, Z, A_H, ΔZ)
        Z_H_min = Int(round(Zₚ) - (NoZperA - 1)/2)
        Z_H_max = Z_H_min + NoZperA - 1
        for Z_H in Z_H_min:Z_H_max
            push!(fragmdomain.A, A_H)
            push!(fragmdomain.Z, Z_H)
            push!(fragmdomain.Value, p_A_Z(Z_H, Zₚ, RMS))
            if A_L != A_H
                Z_L = Z - Z_H
                push!(fragmdomain.A, A_L)
                push!(fragmdomain.Z, Z_L)
                push!(fragmdomain.Value, p_A_Z(Z_L, Z - Zₚ, RMS))
            end
        end
    end
    #Sorting struct in ascending order by mass number
    a = sort(fragmdomain.A)
    z = [fragmdomain.Z[fragmdomain.A .== A] for A in first(a):last(a)]
    z = reduce(vcat, z)
    value = [fragmdomain.Value[fragmdomain.A .== A] for A in first(a):last(a)]
    value = reduce(vcat, value)
    for index_A in eachindex(a)
        fragmdomain.A[index_A] = a[index_A]
        fragmdomain.Z[index_A] = z[index_A]
        fragmdomain.Value[index_A] = value[index_A]
    end
    return fragmdomain
end
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