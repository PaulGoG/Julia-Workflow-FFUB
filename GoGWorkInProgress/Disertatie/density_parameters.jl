#=
Different methods for computing energy-independent density level parameters 
as defined in the Fermi-Gas model:
1. Gilbert-Cameron
2. Egidy-Bucurescu
=#

function density_parameter_Gilbert_Cameron(dGC, A, Z, limInfA_H, limSupA_H)
    a = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in limInfA_H:limSupA_H
        for Z_H in minimum(dGC.n):maximum(dGC.n)
            N_H = A_H - Z_H
            A_L = A - A_H
            Z_L = Z - Z_H
            N_L = A_L - Z_L
            if isassigned(dGC.S_Z[dGC.n .== Z_H], 1) && isassigned(dGC.S_N[dGC.n .== N_H], 1)
                if isassigned(dGC.S_Z[dGC.n .== Z_L], 1) && isassigned(dGC.S_N[dGC.n .== N_L], 1)
                    a_H = A_H * (0.00917*(dGC.S_Z[dGC.n .== Z_H][1] + dGC.S_N[dGC.n .== N_H][1]) + 0.142)
                    a_L = A_L * (0.00917*(dGC.S_Z[dGC.n .== Z_L][1] + dGC.S_N[dGC.n .== N_L][1]) + 0.142)
                    if a_H > 0 && a_L > 0
                        push!(a.x_1, A_H)
                        push!(a.x_2, Z_H)
                        push!(a.y, a_H)
                        push!(a.σ, 0.0)
                        if A_L != A_H
                            push!(a.x_1, A_L)
                            push!(a.x_2, Z_L)
                            push!(a.y, a_L)
                            push!(a.σ, 0.0)
                        end
                    end
                end
            end
        end
    end
    return a
end