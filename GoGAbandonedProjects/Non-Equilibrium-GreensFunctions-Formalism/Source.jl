#Aici scriem functiile de interes

function Lead(P) #Functia ce calculeaza cantitatile la legaturi
    return 𝜮ᴿₚ, 𝜮ⁱₚ, 𝜮ᵒₚ, 𝜞ₚ
end

function Initialise(args) #Initializam 𝜮𝜑 -urile (random?)
    body
    return 𝜮ᴿ𝜑, 𝜮ⁱ𝜑, 𝜮ᵒ𝜑
end

function Conductor(args) #Calculul matricii εI - Hc care se face o singura data, vine la pachet cu calculul ε i.f de moduri si ħvₘ
    body
end

function Green(args) #Calculul functiilor Green
    body
    return Gᴿ, Gᴬ
end

function KineticEq(args) #Rezolvarea ecuatiei cinetice
    body
    return Gⁿ, Gᵖ
end

function Scattering(args) #Recalcularea self-energiilor la imprastiere (CUM?, Discretizarea integralelor)
    body
    return 𝜮ᴿ𝜑, 𝜮ⁱ𝜑, 𝜮ᵒ𝜑
end

function Convergence(args) #Verificarea convergentei (CUM?, cu norma?)
    body
    return Yes or No
end

function Currents(args) #Calculul curentilor (din nou, discretizarea integralelor?)
    body
end
