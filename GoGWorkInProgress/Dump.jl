# Function storage bin
#####

function TKE_A_Z(dy)
    tke = distributie_bidym(Int[], Int[], Float64[], Float64[])
    for A_H in minimum(dy.A_H):maximum(dy.A_H)
        for Z_H in minimum(dy.Z_H[dy.A_H .== A_H]):maximum(dy.Z_H[dy.A_H .== A_H])
            Numarator = 0
            Numitor = sum(dy.Y[(dy.A_H .== A_H) .& (dy.Z_H .== Z_H)])
            Suma_σ² = 0
            # TKE(A, Z) = Σ_(TKE) TKE * Y(A, Z, TKE)/Σ_(TKE) Y(A, Z, TKE)
            # σTKE(A, Z) = [1/Σ_(TKE) Y(A, Z, TKE)] * sqrt[Σ_(TKE) (TKE - TKE(A, Z)) * σY(A, Z, TKE)^2]
            for TKE in minimum(dy.TKE[(dy.A_H .== A_H) .& (dy.Z_H .== Z_H)]):maximum(dy.TKE[(dy.A_H .== A_H) .& (dy.Z_H .== Z_H)])
                Y_A_Z_TKE = sum(dy.Y[(dy.A_H .== A_H) .& (dy.TKE .== TKE) .& (dy.Z_H .== Z_H)])
                Numarator += TKE * Y_A_Z_TKE
            end
            tke_A_Z = Numarator/Numitor
            for TKE in minimum(dy.TKE[(dy.A_H .== A_H) .& (dy.Z_H .== Z_H)]):maximum(dy.TKE[(dy.A_H .== A_H) .& (dy.Z_H .== Z_H)])
                σY_A_TKE = sqrt(sum(dy.σY[(dy.A_H .== A_H) .& (dy.TKE .== TKE) .& (dy.Z_H .== Z_H)].^2))
                Suma_σ² += (TKE - tke_A_Z)^2 * σY_A_TKE^2
            end
            push!(tke.x_1, A_H)
            push!(tke.x_2, Z_H)
            push!(tke.y, tke_A_Z)
            push!(tke.σ, sqrt(Suma_σ²)/Numitor)
        end
    end  
    return tke
end