using

include("Source.jl")

#Raw data
N = #Discretizare
a = #Constanta retelei
ħ = #hbar
m = #masa de banda
#=
VECTORUL LEGATURILOR P ->    Citire CSV sau creere?
Contine valorile fₚ la diferite energii
Pentru fiecare mod m contine energia de cutoff εₘ, nr de unda Kₘ & setul de functii χₘ de
lungime N la fiecare punct de contact dintre conductor si legatura
=#

#Necesitatea definirii potentialului U si a p.mag.vector A discretizat pentru calculul Hc (date de intrare)

#Descrierea interactiilor e-e si e-fonon se face simultan, cuplat, diferit? -> Discretizarea integralelor

t = ħ/(2*m*a)

#Pentru toate legaturile P -> Calcul Lead(P) intr-un for 1->p

#Pentru conductor -> Calcul Conductor()

#Initializam self-energiile de imprastiere cu Initialise()

#=
Deschidem un while in care verificam convergenta
Intr-un for pe toate canalele energetice -> Green()
-> KineticEq() -> Scattering() ->Convergence()

Dupa ce iesim din while calculam curentii -> Currents()
=#
