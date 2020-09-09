using

include("Source.jl")

#Raw data
N = 3 #Discretizare
a = 1 #Constanta retelei

#=
VECTORUL LEGATURILOR ->    Citire CSV sau creere?
Contine valorile fₚ la diferite energii
Pentru fiecare mod m contine energia de cutoff εₘ, nr de unda Kₘ & setul de functii χₘ de
lungime N la fiecare punct de contact dintre conductor si legatura
=#

#Necesitatea definirii potentialului U si a p.mag.vector A discretizat pentru calculul Hc (date de intrare)

#Descrierea interactiilor e-e si e-fonon se face simultan, cuplat, diferit? -> Discretizarea integralelor
