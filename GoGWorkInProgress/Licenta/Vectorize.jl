#=
Aici transformam functiile calculate punctual in vectori
=#

function χ_Q_Durata_Scurta(x, y, z, Pasquill, Suprafata, Tip_Suprafata, t_R)
    χ_Q = zeros(length(x), length(y), length(z))
    for i in 1:length(x)
        for j in 1:length(y)
            for k in 1:length(z)
                χ_Q[i,j,k] = dilutie_instantanee(x[i], y[j], z[k], Pasquill, Suprafata, Tip_Suprafata, t_R)
            end
        end 
    end
    return χ_Q
end
function χ_Durata_Scurta(χ_Q, x, y, z, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
    χ = zeros(length(x), length(y), length(z))
    for i in 1:length(x)
        for j in 1:length(y)
            for k in 1:length(z)
                χ[i,j,k] = χ_scurt(χ_Q[i,j,k], x[i], Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
            end
        end 
    end
    return χ
end
function ω_Scurt(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_R, t_spalare)
    ω = zeros(length(x), length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            ω[i,j] = ω_w_scurt(x[i], y[j], Pasquill, Suprafata, Tip_Aversa, Debit, t_R, t_spalare) + ω_d_scurt(χ_Q[i,j,1], x[i], Pasquill, Suprafata, Tip_Suprafata)
        end 
    end
    return ω
end

function χ_Q_Durata_Prelungita(x, y, Pasquill, Suprafata, Tip_Suprafata)
    χ_Q = zeros(length(x), length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            χ_Q[i,j] = dilutie_durata_prelungita(x[i], y[j], Pasquill, Suprafata, Tip_Suprafata)
        end 
    end
    return χ_Q
end
function χ_Durata_Prelungita(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
    χ = zeros(length(x), length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            χ[i,j] = χ_scurt(χ_Q[i,j], x[i], Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
        end 
    end
    return χ
end
function ω_Prelungit(χ_Q, x, y, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_R, t_spalare)
    ω = zeros(length(x), length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            ω[i,j] = ω_w_scurt(x[i], y[j], Pasquill, Suprafata, Tip_Aversa, Debit, t_R, t_spalare) + ω_d_scurt(χ_Q[i,j], x[i], Pasquill, Suprafata, Tip_Suprafata)
        end 
    end
    return ω
end

function χ_Q_Durata_Lunga(x, y, Suprafata, Tip_Suprafata)
    χ_Q = zeros(length(x), length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            χ_Q[i,j] = dilutie_lunga_durata(x[i], y[j], Suprafata, Tip_Suprafata)
        end 
    end
    return χ_Q
end
function χ_Durata_Lunga(χ_Q, x, y, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
    χ = zeros(length(x), length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            χ[i,j] = χ_lung(χ_Q[i,j], Rotatie(x[i], y[j]), Suprafata, Tip_Suprafata, Tip_Aversa, Debit, Apartenenta_Sector_Cerc(x[i],y[j]), t_spalare)
        end 
    end
    return χ
end
function ω_Lung(χ_Q, x, y, Suprafata, Tip_Suprafata, Tip_Aversa, Debit, t_spalare)
    ω = zeros(length(x), length(y))
    for i in 1:length(x)
        for j in 1:length(y)
            ω[i,j] = ω_w_lung(Rotatie(x[i], y[j]), Suprafata, Tip_Aversa, Debit, t_spalare) + ω_d_lung(χ_Q[i,j], Rotatie(x[i], y[j]), Suprafata, Tip_Suprafata, Apartenenta_Sector_Cerc(x[i], y[j]))
        end
    end
    return ω
end