# Pastram aici un numar de constante vitale pentru calculele programului

const u_10 = 3.5 # Viteza medie vantului la 10 metri deasupra solului in m/s
const Q_0 = 1 # Activitatea poluantului in punctul de emisie in Bq; 1 doar de test!!!
const h  = 150 # Inaltimea cosului de emisie in metri
const h_i = 100 # Inaltimea limitei inferioare a stratului de inversie atmosferica in m; Nu vom implementa strat de inversie!
const n = 16 # Nr de regiuni in care impartim roza vanturilor -> Trebuie sa fie multiplu de 4
const θ_L = 2π/n # Sectorul de unghi din fiecare regiune
const t_R = 2E5 # Timpul de emisie pentru care studiem toate cantitatile in secunde
const Con_Aerodinamic = tan(15*π/180) # limita pozitionarii cladirilor fata de directia vantului
const C = 1 # Intervine la calculul Σ-urilor, normal trebuie sa maximizeze Dilutia
const w_0 = 7 # Viteza de iesire a gazelor pe cos in m/s
const D = 10 # Diametrul interior al cosului in metri
const ρ = 1.2 # Densitatea aerului ambiant
const ρ_0 = 0.7 # Densitatea gazelor emisie; am luat densitatea vaporilor de H2O
const g = 9.8 # Acceleratie gravitationala

const F = (ρ - ρ_0)/ρ * g*w_0*(D/2)^2 # Formula cea mai generala din carte
const F_m = ρ_0/ρ * w_0^2 *(D/2)^2 # Analog F
const C_p = 14 # J(Kg K)
const ΔT_Δz = 9.86E-3 # Gradientul de temperatura in K/m - luat din articol - cuculeanu
const T_0 = 350 # Temperatura gazelor emisie in K
const T = 295 # Temperatura medie a aerului ambiant in K
# const F = (T_0 - T)/T_0 * g*w_0*(D/2)^2 # Alternativa din carte daca seamana gazul emis cu aerul
# const F_m = T/T_0 * w_0^2 *(D/2)^2
const S = g/T * (g/C_p + ΔT_Δz) # Parametru de stabilitate atmosferica

const v_dL_HTO = 0.4E-2 # Minimul vitezei de depunere uscata a HTO
const v_dH_HTO = 0.8E-2 # Maximul vitezei de depunere uscata a HTO
const v_dL_HT = 0.04E-2 # Minimul vitezei de depunere uscata a HT
const v_dH_HT = 0.05E-2 # Maximul vitezei de depunere uscata a HT
const A = 1E-5 # Constanta necesara calculului resuspensiei
const B = 1E-9 # Constanta necesara calculului resuspensiei
const λ_1 = 1E-2 # Constanta necesara calculului resuspensiei
const λ_2 = 2E-5 # Constanta necesara calculului resuspensiei