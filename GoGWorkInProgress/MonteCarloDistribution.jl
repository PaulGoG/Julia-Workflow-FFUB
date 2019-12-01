using Distributions
using Plots

function Distribution(x)
    if x>=0
         return exp(-2*x^2)
    else return exp((-x^2)/2)
    end
end

Np = Int(10e5)
ok=1
contor=1
a = Array{Float64}(undef ,Np)
while ok<=Np
    x = rand(Uniform(-3,3))
    y = rand()
    global contor+=1
    if Distribution(x) > y
        a[ok] = x
        global ok+=1
    end
end

histogram(a, bins=100)

contor/Np
