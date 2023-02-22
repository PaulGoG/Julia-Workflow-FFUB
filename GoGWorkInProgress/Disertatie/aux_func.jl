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

#Function bodies
#Isobaric charge distribution p(Z,A)
function p_A_Z(Z, Z_p, rms_A)
    return 1/(sqrt(2*π) * rms_A) * exp(-(Z - Z_p)^2 /(2* rms_A^2))
end
#
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
#Total Excitation Energy
function Total_excitation_energy(Q, σ_Q, TKE, σ_TKE, Sₙ, σ_Sₙ, Eₙ)
    TXE = Q - TKE + Sₙ + Eₙ
    σ_TXE = sqrt(σ_Q^2 + σ_Sₙ^2 + σ_TKE^2)
    if TXE > 0
        return TXE, σ_TXE
    else 
        return NaN
    end
end
#Construct vectorized fragmentation domain
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
                push!(fragmdomain.Value, p_A_Z(Z_L, Zₚ, RMS))
            end
        end
    end
    return fragmdomain
end
