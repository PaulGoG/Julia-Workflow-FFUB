using Plots
using LaTeXStrings
using LsqFit

# Cod de calcul pentru fitul parametrului Hubble

gr();
cd(@__DIR__); # Adauga calea relativa la folderul de lucru
const c = 2.998 * 10e8;

# Galaxii studiate: NGC 3627. NGC 3368, NGC 3147

function z(λ, λ₀, σ_λ)
    return [(λ - λ₀)/λ₀, sqrt(σ_λ^2 * (λ^2 + λ₀^2)/λ₀^4)]
end
function z_med(z, σ_z)
    Suma_Z = 0.0
    Suma_σ² = 0.0
    for i in eachindex(z)
        Suma_Z += z[i]
        Suma_σ² =+ σ_z[i]^2
    end
    return [Suma_Z/length(z), sqrt(Suma_σ²)/length(z)]
end
function viteza(z_mediu, σ_z)
    return [z_mediu*c, σ_z*c]
end
function distanta(d_et, s_et, s, σ_s_et, σ_s)
    return [d_et*s_et/s, sqrt((d_et*σ_s_et/s)^2 + (d_et*σ_s*s_et/s^2)^2)]
end
function fitare(v, d)
    
end
