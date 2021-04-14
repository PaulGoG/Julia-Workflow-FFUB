#=
Formule diferite pentru calculul dilutiei atmosferice χ/Q
In cazul nostru χ/Q este o matrice NxNxN sau NXN (z=0, la sol)
=#

# Dilutia instantanee sau scurta durata (3 minute -> o ora) -> vector 3D
function dilutie_instantanee(x, y, z, Pasquill)
    # Calculul chestiilor din formula cu functii
    # Se calculeaza pentru o anumita clasa Pasquill
    H = 
    Σ_z = 
    Σ_y = 
    u = 

    χ_Q = 1/(2*π*Σ_y*Σ_z*u) * exp(-y^2/(2*Σ_y^2)) * (exp(-(z-H)^2/(2*Σ_z^2)) + exp(-(z+H)^2/(2*Σ_z^2))) * f_inversie_termica(Σ_z, H, h_i)
    return χ_Q
end

#=
Dilutie durata _prelungita (o ora -> 24 ore)
Depinde de x variabil, y ∈ sector de cerc, z=0 => valori care depind doar de variatia lui x si
Apartententa la un sector unghiular de cerc
Mai exact dependenta explicita va fi doar de x, cea de y fiind folosita pt calculul sectorului unghiular prin arctan
Astfel, valorile o sa fie linii // cu OY marginite de marginile "evantaiului" sectoarelor de cerc
=#
function dilutie_durata_prelungita(x, y)

    # N_sector = floor(atan(abs(y/x)/θ_L)) + 1
    # Se calculeaza pentru o anumita clasa Pasquill
    H = 
    Σ_z = 
    u = 
    χ_Q = (2/π)^0.5 * 1/(Σ_z*u*x*θ_L) * exp(-H^2/(2*Σ_z^2)) * f_inversie_termica(Σ_z, H, h_i)
    return χ_Q
end

# Dilutie lunga_durata (> 24 ore)
# Aici nu mai exista dependenta explicita nici macar de x, totul se considera pe termen foarte lung
# => Sectoare de cerc in care dilutia e cuasi-omogena, depinde de coordonate doar implicit
function dilutie_lunga_durata(x, y)
    # Frecventele se citesc din DataFrame-ul corespunzator
    k = floor(atan(abs(y/x)/θ_L)) + 1
    F_k = first(Freq.F_k[(Freq.Zona_k .== k)])
    F_ki = Freq.F_ki[(Freq.Zona_k .== k)]
    H = 
    Σ_z = 
    u = 
    Parametru = exp(-H^2/(2*Σ_z^2)) * f_inversie_termica(Σ_z, H, h_i)/(Σ_z*u) 
    # Avand in vedere ca toate depind de Pasquill, trebuie ca Parametru sa fie si el vector
    # La fel ca F_ki => Parametru il construim intr-un for alaturi de celalalte calcule pt a trece
    # Prin toate clasele Pasquill cu calculele 
    χ_Q = (2/π)^0.5 * F_k/(x*θ_L) * sum(F_ki .* Parametru)
    return χ_Q
end