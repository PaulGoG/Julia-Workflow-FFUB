#Function bodies used in the main program

#temp structs for testing


struct Distribution{T1 <: Vector{Int}, T2 <: Vector{Float64}}
    A::T1
    Z::T1
    TKE::T2
    Value::T2
    σ::T2
end
