# Pastram aici un numar de constante vitale pentru calculele programului

u_10 = 5 # Viteza medie vantului la 10 metri deasupra solului in m/sa
Q_0 = 1000 # Activitatea poluantului in punctul de emisie in Bq
h  = 70 # Inaltimea cosului de emisie in metri
h_i = 100 # Inaltimea limitei inferioare a stratului de inversie atmosferica in m
n = 16 # Nr de regiuni in care impartim roza vanturilor -> Trebuie sa fie multiplu de 4
θ_L = 2π/n # Sectorul de unghi din fiecare regiune
t_R = 2E5 # Timpul de emisie pentru care studiem toate cantitatile in secunde
C = 1 # Intervine la calculul Σ-urilor, normal trebuie sa maximizeze Dilutia
w_0 = 10 # Viteza de iesire a gazelor pe cos in m/s
D = 20 # Diametrul interior al cosului in metri
ρ = 1.293 # Densitatea aerului ambiant
ρ_0 = 0.269 # Densitatea gazelor emisie (3H)
g = 9.8 # Acceleratie gravitationala
F = 3.7E-5 # Coeficient pe care momentan il scriem asa
C_p = 14 # J(Kg K)
ΔT_z = 9.86E-3 # Gradientul de temperatura in K/m
T_0 = 320 # Temperatura gazelor emisie in K
T = 293 # Temperatura medie a aerului ambiant in K
v_dL_HTO = 0.4E-2 # Minimul vitezei de depunere uscata a HTO
v_dH_HTO = 0.8E-2 # Maximul vitezei de depunere uscata a HTO
v_dL_HT = 0.04E-2 # Minimul vitezei de depunere uscata a HT
v_dH_HT = 0.05E-2 # Maximul vitezei de depunere uscata a HT
A = 1E-5 # Constanta necesara calculului resuspensiei
B = 1E-9 # Constanta necesara calculului resuspensiei
λ_1 = 1E-2 # Constanta necesara calculului resuspensiei
λ_2 = 2E-5 # Constanta necesara calculului resuspensiei