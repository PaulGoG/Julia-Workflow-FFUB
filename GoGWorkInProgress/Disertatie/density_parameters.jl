#=
Different methods for computing energy-independent density level parameters 
as defined in the Fermi-Gas model:
1. Gilbert-Cameron
2. Egidy-Bucurescu
=#
#####
function density_parameter_Gilbert_Cameron(A, Z, dGC)
    N = A - Z
    if isassigned(dGC.S_Z[dGC.n .== Z], 1) && isassigned(dGC.S_N[dGC.n .== N], 1)
        a = A * (0.00917*(dGC.S_Z[dGC.n .== Z][1] + dGC.S_N[dGC.n .== N][1]) + 0.142)
        if a > 0
            return a
        else
            return NaN
        end
    else 
        return NaN
    end
end
function density_parameter_Egidy_Bucurescu(A, Z, dm)
    if isassigned(dm.D[(dm.A .== A) .& (dm.Z .== Z)], 1) && isassigned(dm.D[(dm.A .== A+2) .& (dm.Z .== Z+1)], 1) && isassigned(dm.D[(dm.A .== A-2) .& (dm.Z .== Z-1)], 1)
        D = dm.D[(dm.A .== A) .& (dm.Z .== Z)][1]*1e-3
        W_exp = Z*Dᵖ + (A-Z)*Dⁿ - D
        a_sim = A*(a_sim_1 - a_sim_2 *A^(-1/3))
        η = (A - 2*Z)/A
        W_LDM = a_v*A - a_s*A^(2/3) - a_c*Z^2 *A^(-1/3) - a_sim*η^2
        δW₀ = W_LDM - W_exp
        D_plus = dm.D[(dm.A .== A+2) .& (dm.Z .== Z+1)][1]*1e-3
        D_minus = dm.D[(dm.A .== A-2) .& (dm.Z .== Z-1)][1]*1e-3
        P_d = (D_plus - 2*D + D_minus)/2
        δW = δW₀ + P_d/2
        a = (p_eb_1 + p_eb_2*δW) *A^p_eb_3
        if a > 0
            return a
        else 
            return NaN
        end
    else
        return NaN
    end
end
function density_parameter(density_parameter_type, A, Z, density_parameter_datafile)
    if density_parameter_type == "GC"
        a = density_parameter_Gilbert_Cameron(A, Z, density_parameter_datafile)
    elseif density_parameter_type == "BSFG"
        a = density_parameter_Egidy_Bucurescu(A, Z, density_parameter_datafile)
    end
    return a
end