# Numerical calculation of a simple integral from experimental data using trapezoid method

using Trapz
using Plots
using Interpolations

y = [0,0.0953,0.182,0.262,0.336,0.405,0.470,0.531,0.588,0.642,0.693]
y = [1,1.067,1.13,1.232,1.356,1.485,1.624,1.781,1.925,2.22,2.275]
x = collect(1:0.1:2)
I = trapz(x,y)

itp = interpolate(y , BSpline(Cubic(Reflect(OnGrid())))) # Domeniu arbitrar
sitp = scale(itp, 1:0.1:2) #Rescalare
plot(itp)

a = collect(1:0.01:2)
plot(a, sitp.(a))


b = itp.(a)
scatter(a,b)
scatter!(x,y)
