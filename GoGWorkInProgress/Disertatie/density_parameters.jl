#=
Different methods for computing energy-independent density level parameters 
as defined in the Fermi-Gas model:
1. Gilbert-Cameron
2. Egidy-Bucurescu
=#

function density_parameter_Gilbert_Cameron(A, Z, A_H, Z_H, dGC)
    N_H = A_H - Z_H
    A_L = A - A_H
    Z_L = Z - Z_H
    N_L = A_L - Z_L
    if isassigned(dGC.S_Z[dGC.n .== Z_H], 1) && isassigned(dGC.S_N[dGC.n .== N_H], 1)
        if isassigned(dGC.S_Z[dGC.n .== Z_L], 1) && isassigned(dGC.S_N[dGC.n .== N_L], 1)
            a_H = A_H * (0.00917*(dGC.S_Z[dGC.n .== Z_H][1] + dGC.S_N[dGC.n .== N_H][1]) + 0.142)
            a_L = A_L * (0.00917*(dGC.S_Z[dGC.n .== Z_L][1] + dGC.S_N[dGC.n .== N_L][1]) + 0.142)
            if a_H > 0 && a_L > 0
                return a_L, a_H
            else
                return NaN
            end
        else 
            return NaN
        end
    else 
        return NaN
    end
end