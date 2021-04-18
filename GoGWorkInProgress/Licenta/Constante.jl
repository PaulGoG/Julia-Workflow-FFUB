# Pastram aici un numar de constante vitale pentru calculele programului

const u_10 = 3.5 # Viteza medie vantului la 10 metri deasupra solului in m/s
const Q_0 = 1E16 # Activitatea poluantului in punctul de emisie in Bq
const h  = 150 # Inaltimea cosului de emisie in metri
const n = 16 # Nr de regiuni in care impartim roza vanturilor -> Trebuie sa fie multiplu de 4!
const θ_L = 2*π/n # Sectorul de unghi din fiecare regiune
const Sector_Cerc = tan(θ_L / 2)
const C = 2 # Intervine la calculul Σ-urilor, trebuie sa maximizeze dilutia

const w_0 = 5 # Viteza de iesire a gazelor pe cos in m/s
const D = 15 # Diametrul interior al cosului in metri
const ρ = 1.225 # Densitatea aerului ambiant
const ρ_0 = 0.85 # Densitatea gazelor emisie; am luat densitatea vaporilor de H2O
const g = 9.8 # Acceleratie gravitationala

const F = (ρ - ρ_0)/ρ * g*w_0*(D/2)^2 # Formula cea mai generala din carte
const F_m = ρ_0/ρ * w_0^2 *(D/2)^2 # Analog F
const C_p = 1880 # J(Kg K)
const ΔT_Δz = 9.86e-3 # Gradientul de temperatura in K/m - luat din articol - cuculeanu
const T_0 = 350 # Temperatura gazelor emisie in K
const T = 295 # Temperatura medie a aerului ambiant in K
# const F = (T_0 - T)/T_0 * g*w_0*(D/2)^2 # Alternativa din carte daca seamana gazul emis cu aerul
# const F_m = T/T_0 * w_0^2 *(D/2)^2
const S = g/T * (g/C_p + ΔT_Δz) # Parametru de stabilitate atmosferica

const λ_i = 1.784E-9 # Constanta de dezintegrare a tritiului in 1/s
const v_dL_HTO = 0.4e-2 # Minimul vitezei de depunere uscata a HTO
const v_dH_HTO = 0.8e-2 # Maximul vitezei de depunere uscata a HTO
const v_dL_HT = 0.04e-2 # Minimul vitezei de depunere uscata a HT
const v_dH_HT = 0.05e-2 # Maximul vitezei de depunere uscata a HT

const A = 1E-5 # Constanta necesara calculului resuspensiei in m^-1
const B = 1E-9 # Constanta necesara calculului resuspensiei in m^-1
const λ_1 = 1E-2 # Constanta necesara calculului resuspensiei in zile^-1
const λ_2 = 2E-5 # Constanta necesara calculului resuspensiei in zile^-1