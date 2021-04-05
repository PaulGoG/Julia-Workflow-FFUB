# Function bodies that are going to be used in the calculations

# Ma gandesc ca aici sus sa fie citite Data Frames -> datele tabelate pe care sa le acceseze functiile


# Functii corectie inaltime efectiva → 2 corectii simple & directe
# A treia corectie in cazul portantei re-apeleaza functia de calcul a dispersiei care la randul sau se calculeaza cu H
# Nu mai folosim momentan strat de inversie!!! => scapam de recursivitate
# Nu inteleg diferenta intre faza tranzitie vs faza finala


# Functii calcul dispersii 
function σ_z(x, Pasquill)
    #inserare tabel 1 si tabel 2 ca dependinte in functii deja citite!. Clasa Pasquill A->F cu indicii 1->6

    g = (a_1*x^b_1)/(1+a_2*x^b_2)

    if z_0 > 0.1
        F = log(c_1 * x^d_1 * (1+(c_2*x^d_2^(-1))))
    else
        F = log((c_1*x^d_1)/(1+c_2*x^d_2))
    end

    return g*F
end

function σ_y(x, Pasquill, timp)
    #Tabel pt c_3 tabel 3
    σ_y = (c_3*x)/(1+0.0001*x)^0.5
    if (t < 600)
         return σ_y
    else
         return σ_y * (t/600)^0.2
    end

end

# Functie corectii dispersii (entangled cu corectia de portanta)
function Σ_y(x,y?)
    #Necesita structura cladirilor din zona in ceva matricea
    #Calcul dsursa-cladire = sqrt(x^2+y^2) [vector]
    #Pastram doar valorile pt care ds-c < 3*zcladire & arctg(modul(y)/modul(x)<5 grade sa fie apropiata cladirea de OX)
    #Dupa triaj calculam 2 marimi mediate pentru formule 
    #H_cladire = suma(z/dsc)/suma(1/dsc) & analog Acladire
    #Ne trebuie neaparat H

    if H >= 2.5*H_cladire 
        Σ_y = σ_y
    else
        Σ_max = (σ_y^2 + C*A_cladire/π)^0.5
        if H < H_cladire 
            Σ_y = Σ_max
        else 
            Σ_y = Σ_max - (H-H_cladire)/(1.5 * H_cladire) * (Σ_max - σ_y)
        end
    end
    return Σ_y
end

function Σ_z(?)
    #Analog exact dupa Σ_y
    #Ma gandesc ca partea cu evaluarea d-urilor la cladiri sa fie facuta separat, o singura data in alta functie...
    #Sau chiar sa fie rulata la inceputul blocului astuia de functii, sau functie in main, vedem...
end

# Functie calcul strat de inversie termica (aici am nevoie de detalii si indrumare teoretice)
# Inversia ∃ pe tot spatiul?, sau e data de masurari meteo?
# Necesita H
# Momentan rulam tot cu f = 1 si vedem dupa aceea daca mai complicam...

# Calcul factor DEC fara descendent; T se dezintegreaza in He stabil si stim lambda din literatura
function DEC(x)
    return exp(-λ*x/u)
end

# Calcul factor DEP -> DOAR DEPUNERE USCATA!, depunerea umeda necesita masuratori meteo precise
function DEP_scurt()
    #Integrarea numerica
end

function DEP_lung()
    #α si dupa integrarea numerica
end
# Calcul concentratie integrata in timp χ
function χ()
    #χ_Q * DEC * DEP * Q_0 sau Q_med
end

# Calcul depuneri uscate
function ω_d()
    #χ * v_d care e in carte pt HTO sau HT
end

# Calcul Resuspensii
function Resuspensie()
    K = A*exp(-λ_1*t) + B*exp(-λ_2*t)
    #Valori de referinta in concentratie
    #Returneaza K-ul si mai apoi cu el se calculeaza punctual 
    #Resuspensia ca fiind K*ω => reprezentare grafica 
end