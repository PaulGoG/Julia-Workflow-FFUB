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
        if d <= 3*Cladiri.z[i] && atan(abs(Cladiri.y[i]/Cladiri.x[i])) < 20*π/180
            d_sc[i] = d
            H_cl += Cladiri.z[i]/d
            A_cl += Cladiri.Arie_transversala[i]/d
        end        
    end
    filter!(x -> x!=0, d_sc)
    Suma = sum(1 ./d_sc)
    return H_cl/Suma, A_cl/Suma
end