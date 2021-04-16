#=
Formule diferite pentru calculul dilutiei atmosferice χ/Q
In cazul nostru χ/Q este o matrice NxNxN sau NXN (z=0, la sol)
!Simulam pana de poluant pentru care x>0 pe directia vantului!
=#

#=
Dilutia instantanee sau scurta durata (3 minute -> o ora) -> vector 3D
A nu se apela in punctul de emisie sau impotriva vantului (x<=0)!
=#
function dilutie_instantanee(x, y, z, Pasquill, Suprafata, Tip_Suprafata)
    if x > 0
    H = H_final(x, Pasquill, Suprafata)
    Σz = Σ_z(x, Pasquill, Suprafata, Tip_Suprafata)
    Σy = Σ_y(x, Pasquill, Suprafata)
    u = u_z(H, Pasquill, Suprafata)
    return 1/(2*π*Σy*Σz*u) * exp(-y^2/(2*Σy^2)) * (exp(-(z-H)^2/(2*Σz^2)) + exp(-(z+H)^2/(2*Σz^2)))
    end
    return 0.0
end

#=
Dilutie durata _prelungita (o ora -> 24 ore)
Depinde de x variabil, y ∈ sector de cerc, z=0 => valori care depind doar de variatia lui x si
Apartententa la un sector unghiular de cerc
Mai exact dependenta explicita va fi doar de x, cea de y fiind folosita pt calculul sectorului unghiular prin arctan
Astfel, valorile o sa fie linii // cu OY marginite de marginile "evantaiului" sectorului de cerc
=#
function dilutie_durata_prelungita(x, y, Pasquill, Suprafata, Tip_Suprafata)
    if x > 0 && abs(y/x) <= Sector_Cerc
    H = H_final(x, Pasquill, Suprafata)
    Σz = Σ_z(x, Pasquill, Suprafata, Tip_Suprafata)
    u = u_z(H, Pasquill, Suprafata)
    return (2/π)^0.5 * 1/(Σz*u*x*θ_L) * exp(-H^2/(2*Σz^2))
    end
    return 0.0
end

# Dilutie lunga_durata (> 24 ore)
# Aici nu mai exista dependenta explicita nici macar de x, totul se considera pe termen foarte lung
# => Sectoare de cerc in care dilutia e cuasi-omogena, depinde de coordonate doar implicit
function dilutie_lunga_durata(x, y, Suprafata, Tip_Suprafata)
        k = Apartenenta_Sector_Cerc(x, y)
        xrotit = Rotatie(x,y)
        F_k = Freq.F_k[(Freq.Zona_k .== k)][1]
        Suma = 0.0
            for i in 1:6
                Pasquill = IntegertoString(i)
                H = H_final(xrotit, Pasquill, Suprafata)
                Σz = Σ_z(xrotit, Pasquill, Suprafata, Tip_Suprafata)
                u = u_z(H, Pasquill, Suprafata)
                F_ki = Freq.F_ki[(Freq.Zona_k .== k) .& (Freq.Clasa_Pasquill .== Pasquill)][1]
                Suma = Suma + F_ki * exp(-H^2/(2*Σz^2))/(Σz*u) 
            end
        return (2/π)^0.5 * F_k * Suma/(xrotit*θ_L)
end