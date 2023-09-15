using QuadGK
using Plots

R = 30e-3
P_l = 5.1e-8
N_c = 64e3
g = 0.23
y = 15
x_0 = 50.0
x_range = collect(-x_0:1.0:x_0)

function Caz_1(x)
    return (abs(x)/sqrt(x^2 + y^2)) * exp(-g * sqrt(x^2 + y^2))
end

function Caz_2(x)
    if x <= 0
        return (abs(x)^2/(x^2 + y^2)) * exp(-g * sqrt(x^2 + y^2))
    else
        return 0
    end
end

I_1 = g^2/(2*π) * first(quadgk(Caz_1, -x_0, x_0))
I_2 = g^2/2 * first(quadgk(Caz_2, -x_0, x_0))

P_1 = P_l * N_c * π * R^2 * I_1;
P_2 = P_l * N_c * π * R^2 * I_2;

plt = plot(x_range, Caz_2);
#plt = plot!(x_range, Caz_2);
display(plt)

println(P_2)