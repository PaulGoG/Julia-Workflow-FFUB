#=
Aici transformam functiile calculate punctual in vectori
=#

function χ_Q_χ_Durata_Scurta(x, y, z, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
    χ_Q = zeros(length(x), length(y), length(z))
    χ = zeros(length(x), length(y), length(z))
    for i in 1:length(x)
        for j in 1:length(y)
            for k in 1:length(z)
                χ_Q[i,j,k] = dilutie_instantanee(x[i],y[j],z[k], Pasquill, Suprafata, Tip_Suprafata)
                χ[i,j,k] = χ_scurt(χ_Q[i,j,k], x[i], Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
            end
        end 
    end
    return χ_Q, χ
end
function ω_K_Scurt_Prelungit(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)
    ω = zeros(length(x), length(y))
    K = zeros(length(x), length(y))
    for i in 1:length(x)
        K[i,:] .= Resuspensie_scurt(x[i], Pasquill, Suprafata)
        for j in 1:length(y)
            ω[i,j] = ω_w_scurt(x[i], y[j], Pasquill, Suprafata, Tip_Aversa, Debit) + ω_d_scurt(χ_Q[i,j,1], x[i], Pasquill, Suprafata, Tip_Suprafata)
        end 
    end
    return χ_Q, χ, ω, K
end