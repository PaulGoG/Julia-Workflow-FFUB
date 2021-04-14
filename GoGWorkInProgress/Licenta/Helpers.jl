# Corpul functiilor principale din program

# Conversie String - Integer si vice-versa pentru lucrul cu clase Pasquill
function StringtoInteger(String)
    if String == "A"
        return 1
    elseif String == "B"
        return 2
    elseif String == "C"
        return 3
    elseif String == "D"
        return 4
    elseif String == "E"
        return 5
    elseif String == "F"
        return 6
    end
    return nothing
end
function IntegertoString(Integer)
    if Integer == 1
        return "A"
    elseif Integer == 2
        return "B"
    elseif Integer == 3
        return "C"
    elseif Integer == 4
        return "D"
    elseif Integer == 5
        return "E"
    elseif Integer == 6
        return "F"
    end
    return nothing
end

#=
Calculul vitezei vantului pe directia OX la o anumita inaltime z
Suprafata poate fi "Apa", "Agricol", "Padure_Oras"
Pasquill poate fi orice string de la "A" la "F"
Valoarea este constanta dupa inaltimea de 200 m
=#

# Calculul u(z) intr-un punct -> Depinde de Pasquill
function u_z(z, Pasquill, Suprafata) 
    m = T_4.m[(T_4.Clasa_Pasquill .== Pasquill) .& (T_4.Tip_Suprafata .== Suprafata)]
    if z <= 200
        return u_10*(z/10.0)^m[1]
    end
    return u_10*(200/10.0)^m[1]
end
# Constructia vectorului de valori u(z)
function Construct_Vector_u_z(z, Pasquill, Suprafata)
    A = zeros(length(z))
    for i in 1:length(z)
        A[i] = u_z(z[i], Pasquill, Suprafata)
    end
    return A
end
# Calculul unui u mediat pe clasele Pasquill
function u_mediu_z(z, Suprafata)
    u = 0.0
    for i in 1:6
        u = u + u_z(z, IntegertoString(i), Suprafata)
    end
    return u/6
end
# Constructia vectorului de valori u_med
function Construct_Vector_u_mediu_z(z, Suprafata)
    A = zeros(length(z))
    for i in 1:length(z)
        A[i] = u_mediu_z(z[i], Suprafata)
    end
    return A
end

# Functii corectie inaltime efectiva → curenti de aer descendenti & antrenare in cavitatea aerodinamica a cladirilor
function H_1(Suprafata)
    u = u_mediu_z(h, Suprafata)
    if w_0 < 1.5*u
        return h - 2*(1.5-w_0/u)*D
    else
        return h
    end
end
function H_2(Suprafata)
    H1 = H_1(Suprafata)
    H_cladire = Echivalent_Cladire()[1]
    if H1 < H_cladire 
        return 0
    else
        if H1 > 2.5 * H_cladire 
            return H1
        else 
            if u_mediu_z(h, Suprafata) < 5
                return H1
            else
                return H1 - (1.5*H_cladire - 0.6*H1)
            end
        end
    end
end

#=
Implementarea cladirilor -> conform teoriei furnizam un H_cladire si A_cladire care reflecta
Contributia fiecarei cladiri normata la apropierea ei de cos (inversul distantei)
Se considera oricum doar cladirile pentru care x^2+y^2 <= (3*h_cladire)^2 &
Din cadranele I si IV trigonometrice verificate cu arctan
=#

function Echivalent_Cladire()
    d_sc = zeros(length(Cladiri[:,1]))
    H_cl = 0.0
    A_cl = 0.0
    for i in 1:length(Cladiri[:,1])
        d = sqrt(Cladiri.x[i]^2 + Cladiri.y[i]^2)
        if d <= 3*Cladiri.z[i] && Cladiri.y[i]/Cladiri.x[i] <= Con_Aerodinamic && Cladiri.y[i]/Cladiri.x[i] >= 0 
            d_sc[i] = d
            H_cl += Cladiri.z[i]/d
            A_cl += Cladiri.Arie_transversala[i]/d
        end        
    end
    filter!(x -> x != 0, d_sc)
    Suma = sum(1 ./d_sc)
    return H_cl/Suma, A_cl/Suma
end

#=
Implementarea corectiilor suprainaltarii penei de poluant din cauza
Impulsului gazelor la iesirea din cos & a portantei gazelor din cauza
Diferentei de temperatura fata de mediul ambiant
=#

function X_0()
    if F < 55
        return 14*F^(5/8)
    end
        return 34*F^(2/5)
end
const x_0 = X_0() # Parametru de care avem nevoie la corectia Δh_b

# Corectia efectului de portanta
function Δh_b_final(Pasquill, Suprafata)
    u = u_z(H_2(Suprafata), Pasquill, Suprafata)
    dh_b_stabil = 2.6*(F/(u*S))^(1/3)
    dh_b_neutru = 1.6*F^(1/3) * (3.5*x_0)^(2/3) /u
    dh_b_calm = 5.0 * F^(1/4) * S^(-3/8)
    return min(dh_b_stabil, dh_b_neutru, dh_b_calm)
end
function Δh_b(x, Pasquill, Suprafata)
    u = u_z(H_2(Suprafata), Pasquill, Suprafata)
    hbfinal = Δh_b_final(Pasquill, Suprafata)
    hbtranzitie = 1.6*F^(1/3) * x^(2/3) / u
    if x < 3.5*x_0 && hbtranzitie <= hbfinal
        return hbtranzitie
    end
    return hbfinal
end

# Corectia efectului de impuls
function Δh_m_final(Pasquill, Suprafata)
    u = u_z(H_2(Suprafata), Pasquill, Suprafata)
    dh_m_neutru = 1.5*w_0*D/u
    dh_m_calm = 4*(F_m/S)^(1/4)
    dh_m_stabil = 1.5 * (F_m/u)^(1/3) * S^(-1/6)
    return min(dh_m_neutru, dh_m_calm, dh_m_stabil)
end
function Δh_m(x, Pasquill, Suprafata)
    u = u_z(H_2(Suprafata), Pasquill, Suprafata)
    hmfinal = Δh_m_final(Pasquill, Suprafata)
    hmtranzitie = 1.89*(w_0^2 * D/(u*(w_0 + 3*u)))^(2/3) * x^(1/3)
    if hmtranzitie <= hmfinal
        return hmtranzitie
    end
    return hmfinal
end

# Corectia combinata cu formula semiempirica
function Δh_mb(x, Pasquill, Suprafata)
    u = u_z(H_2(Suprafata), Pasquill, Suprafata)
    hmbfinal = Δh_m_final(Pasquill, Suprafata) + Δh_b_final(Pasquill, Suprafata)
    hmbtranzitie = 3^(1/3) * (F_m*x/((1/3 + u/w_0)^2 * u^2) + F*x^2/(0.5 * u^3))^(1/3)
    if hmbtranzitie <= hmbfinal
        return hmbtranzitie
    end
    return hmbfinal
end

# Valoarea finala-corectata a inaltimii efective de emisie
function H_final(x, Pasquill, Suprafata)
    hbfinal = Δh_b_final(Pasquill, Suprafata)
    hmfinal = Δh_m_final(Pasquill, Suprafata)
    if abs(hbfinal - hmfinal)*2/(hbfinal + hmfinal) <= 0.1
        return  H_2(Suprafata) + Δh_mb(x, Pasquill, Suprafata)       
    else
        if hmfinal > hbfinal
            return H_2(Suprafata) + Δh_m(x, Pasquill, Suprafata)
        else
            return H_2(Suprafata) + Δh_b(x, Pasquill, Suprafata)
        end
    end
end


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
function Σ_y(x,y)
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

function Σ_z()
    #Analog exact dupa Σ_y
    #Ma gandesc ca partea cu evaluarea d-urilor la cladiri sa fie facuta separat, o singura data in alta functie...
    #Sau chiar sa fie rulata la inceputul blocului astuia de functii, sau functie in main, vedem...
end

# Functie calcul strat de inversie termica (aici am nevoie de detalii si indrumare teoretice)
# Inversia ∃ pe tot spatiul?, sau e data de masurari meteo?
# Necesita H
# Momentan rulam tot cu f = 1 si vedem dupa aceea daca mai complicam...
function f_inversie_termica(Σ_z, H, h_i)
    return 1 
end

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
    t_zile = (t_R - x/u) * 86400
    K = A*exp(-λ_1*t_zile) + B*exp(-λ_2*t_zile) 
    #Valori de referinta in concentratie
    #Returneaza K-ul si mai apoi cu el se calculeaza punctual 
    #Resuspensia ca fiind K*ω => reprezentare grafica 
end