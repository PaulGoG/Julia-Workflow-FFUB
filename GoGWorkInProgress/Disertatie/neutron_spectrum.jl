#Iteratively compute total neutron spectrum in Laboratory frame from main DSE data

struct Neutron_spectrum{T <: Vector{Float64}} <: AbstractDistribution
    ε::T
    Value::T
end

function Neutron_spectrum(ε::Float64, T::Float64)
    return ε*exp(-ε/T)/T^2
end
function Neutron_spectrum(ε::Float64, α::Float64, T::Float64)
    (ε + α*sqrt(ε))*exp(-ε/T) /((sqrt(T) +α*sqrt(π)/2) *T^(3/2))
end