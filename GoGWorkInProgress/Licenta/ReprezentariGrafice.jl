# Singura chestie care va varia in timp in mod explicit va fi resuspensia din cauza lui K, restul sunt
# valori statice/mediate -> heatmap sau contourf sau surface

function DurataScurta(dim_transversal, dim_vertical, Pasquill, Suprafata, Tip_Suprafata, Tip_Aversa, Debit)

    x = collect(1:1:dim_transversal)
    y = collect(-(dim_transversal/2):1:(dim_transversal/2))
    z = collect(0:1:dim_vertical)

    χ_Q = zeros(dim_transversal + 1, dim_transversal + 1, dim_vertical + 1)
    for i in 1:(dim_transversal + 1)
        for j in 1:(dim_transversal + 1)
            for k in 1:(dim_vertical + 1)
                χ_Q[i,j,k] = dilutie_instantanee(x[i],y[j],z[k], Pasquill, Suprafata, Tip_Suprafata)
            end
        end 
    end

    ω = zeros(dim_transversal + 1, dim_transversal + 1)
    for i in 1:length(x)
        ω[i,:] = ω_w_scurt(x[i], Pasquill, Suprafata, Tip_Aversa, Debit)
        for j in 1:length(y)
            for k in 1:length(z)
                
            end
        end 
    end

end