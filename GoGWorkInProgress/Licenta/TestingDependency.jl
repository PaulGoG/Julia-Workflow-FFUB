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
end
function u_z(z, Pasquill, Suprafata)
    m = T_4.m[(T_4.Clasa_Pasquill .== Pasquill) .& (T_4.Tip_Suprafata .== Suprafata)]
    if z <= 200
        return u_10*(z/10.0)^m[1]
    end
    return u_10*(200/10.0)^m[1]
end
function Construct_Vector_u_z(z, Pasquill, Suprafata)
    A = zeros(length(z))
    for i in 1:length(z)
        A[i] = u_z(z[i], Pasquill, Suprafata)
    end
    return A
end
function u_mediu_z(z, Suprafata)
    u = 0.0
    for i in 1:6
        u = u + u_z(z, IntegertoString(i), Suprafata)
    end
    return u/6
end
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

function X_0()
    if F < 55
        return 14*F^(5/8)
    end
        return 34*F^(2/5)
end
const x_0 = X_0() # Parametru de care avem nevoie la corectia Δh_b


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

function Δh_mb(x, Pasquill, Suprafata)
    u = u_z(H_2(Suprafata), Pasquill, Suprafata)
    hmbfinal = Δh_m_final(Pasquill, Suprafata) + Δh_b_final(Pasquill, Suprafata)
    hmbtranzitie = 3^(1/3) * (F_m*x/((1/3 + u/w_0)^2 * u^2) + F*x^2/(0.5 * u^3))^(1/3)
    if hmbtranzitie <= hmbfinal
        return hmbtranzitie
    end
    return hmbfinal
end

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