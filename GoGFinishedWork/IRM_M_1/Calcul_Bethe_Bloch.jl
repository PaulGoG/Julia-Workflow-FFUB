termen_constant = 2*pi*6.023e23* (2.817e-15)^2 *0.511;

termen_tinta = 2.33e6 *14/28;

E_cin_alfa = 5;
M_alfa_C_patrat = 4.0015065 * 931.5;
gamma = E_cin_alfa/M_alfa_C_patrat + 1;
beta = sqrt(1 - 1/gamma^2);
termen_proiectil = 2^2 /beta^2;

potential_interactie = (9.76*14 + 58.8*14^(-0.19)) * 1e-6;
termen_paranteza = log(2*0.511*gamma^2 *beta^2 *2*0.511*gamma^2 *beta^2/potential_interactie^2) -2*beta^2;

BB = termen_constant*termen_tinta*termen_proiectil*termen_paranteza; # In MeV/m
println("dE/dx = ", BB/100, " MeV/cm")