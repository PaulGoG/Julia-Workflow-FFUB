# Pastram aici un numar de constante vitale pentru calculele programului

const u_10 = 5 # Viteza medie vantului la 10 metri deasupra solului in m/s
const Q_0 = 1000 # Activitatea poluantului in punctul de emisie in Bq
const h  = 150 # Inaltimea cosului de emisie in metri
const h_i = 100 # Inaltimea limitei inferioare a stratului de inversie atmosferica in m
const n = 16 # Nr de regiuni in care impartim roza vanturilor -> Trebuie sa fie multiplu de 4
const θ_L = 2π/n # Sectorul de unghi din fiecare regiune
const t_R = 2E5 # Timpul de emisie pentru care studiem toate cantitatile in secunde
const C = 1 # Intervine la calculul Σ-urilor, normal trebuie sa maximizeze Dilutia
const w_0 = 7 # Viteza de iesire a gazelor pe cos in m/s
const D = 20 # Diametrul interior al cosului in metri
const ρ = 1.293 # Densitatea aerului ambiant
const ρ_0 = 0.269 # Densitatea gazelor emisie (3H)
const g = 9.8 # Acceleratie gravitationala
const F = 3.7E-5 # Coeficient pe care momentan il scriem asa
const C_p = 14 # J(Kg K)
const ΔT_Δz = 9.86E-3 # Gradientul de temperatura in K/m
const T_0 = 320 # Temperatura gazelor emisie in K
const T = 293 # Temperatura medie a aerului ambiant in K
const v_dL_HTO = 0.4E-2 # Minimul vitezei de depunere uscata a HTO
const v_dH_HTO = 0.8E-2 # Maximul vitezei de depunere uscata a HTO
const v_dL_HT = 0.04E-2 # Minimul vitezei de depunere uscata a HT
const v_dH_HT = 0.05E-2 # Maximul vitezei de depunere uscata a HT
const A = 1E-5 # Constanta necesara calculului resuspensiei
const B = 1E-9 # Constanta necesara calculului resuspensiei
const λ_1 = 1E-2 # Constanta necesara calculului resuspensiei
const λ_2 = 2E-5 # Constanta necesara calculului resuspensiei