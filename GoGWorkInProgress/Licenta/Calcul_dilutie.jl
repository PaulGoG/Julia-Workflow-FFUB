# Different formulae for calculating the dillution χ/Q
# In our case χ/Q is an NxNxN matrix representing the dillution in the whole Euclidian space

# Dilutia instantanee -> vector tri-dimensional
function dilutie_instantanee(x,y,z)
    χ_Q = 1/(2*π*Σ_y*Σ_z*u) * exp(-y^2/(2*Σ_y^2)) * (exp(-(z-H)^2/(2*Σ_z^2)) + exp(-(z+H)^2/(2*Σ_z^2))) * f_inversie_termica(Σ_z, H, h_i)
end

# Dilutie d_p -> x variabil, y ∈ sector de cerc, z=0 => valori care depind doar de variatia lui x si apartententa la un sector unghiular de cerc
function dilutie_durata_prelungita(x, y)

    θ_L = 2*π/16
    N_sector = modf(atan(y, x)/θ_L)[2]

    χ_Q = (2/π)^0.5 * 1/(Σ_z*u*x*θ_L) * exp(-H^2/(2*Σ_z^2)) * f_inversie_termica(Σ_z, H, h_i)
    
end

# Dilutie l_d
function dilutie_lunga_durata(x, y)

    θ_L = 2*π/16
    k = modf(atan(y, x)/θ_L)[2]

    χ_Q = (2/π)^0.5 * 1/(Σ_z*u*x*θ_L) * exp(-H^2/(2*Σ_z^2)) * f_inversie_termica(Σ_z, H, h_i)
    
end