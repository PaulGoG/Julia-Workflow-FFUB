using DifferentialEquations
using Plots
using CSV
using LsqFit

const G = 6.67408*1e-11        #Cosmological constant

"""
        Friedmann(du, u, p, t)

    Creates the Friedmann equations system of differential equations:
    dH/dt = -H^2 -𝜺(1+3w) +Λ
    da/dt = Ha sau dz/dt = -H(1+z)
    d𝜺/dt = -3H𝜺(1+w) +Λ
    With the constraint 𝜺a^3(1+w) = constant and non-singular initial values

"""
function FriedmannTimeDependentNeutralFluid(du, u, p, t)
    du[1] = -u[1]^2 -(4*π*G/3)*u[3]*(1+3*p[1]) +p[2]
    du[2] = u[1]*u[2]
    du[3] = -3*u[1]*u[3]*(1+p[1]) +p[2]
end

"""
    dH/dz = 1/(z+1) [H + 4πG/3H (𝜺 + 3p)]
    da/dz = - a/(z+1)
    d𝜺/dz = 3/(z+1) (𝜺 + p)

    Cu cazurile pentru p = f(𝜺):
        I Quintesence --> p = w𝜺 ; w∈(-1,-1/3)
        II Chapligyn gas --> p = -A/𝜺^α ;  A>0 & α∈(0,1)
        III Modified Chapligyn gas --> p = w𝜺 -A/𝜺^α cu constrangerile de mai sus
"""


# Plot pentru H in functie de z analitic (constante 1)
x = rand(0.1:2., 1, 50)
y = rand(1,50)
for i =1:length(x)
    y[i]=  sqrt(fit.param[1]*x[i] + fit.param[2]*x[i]^5)
end
plot(x, y, legend=:false)




function HDependentOnZ(du, u, p, t)
    du[1] = u[1]/(t+1) + (4*π*G*(1+3*p[1])* p[2]/3) * (1/u[1]) * (t+1)^(3*(1 +p[1]) -1) #p[1] e w, p[2] e 𝜺0
end

f = HDependentOnZ
tspan = (0.07,2.)
u0 = [68.]
p0 = [1/3, 1e12]
prob = ODEProblem(f, u0, tspan, p0)
sol = solve(prob)
plot(sol)
str = Array(sol)

# Rezolvarea numerica pentru H in functie de z + plotare
# Ecuatia lui H ca functie de z: dH/dz = H/(z+1) + 4πG(1+3w)𝜺0/3 * 1/H * (z+1)^[3(1+w)-1]

h(H, w, z) =  H/(z+1) + (4*π*G*(1+3*w)*1e12/3) * (1/H) * (z+1)^(3*(1+w)-1)
H0 = 68.
w = 1/3
zspan = (0.07,2.)
prob = ODEProblem(h, H0, zspan, w)
sol = solve(prob)
plot(sol)
x = sol.t
y = sol.u
plot(x,y)





#  Rezolvarea numerica a sistemului de ecuatii Friedmann ce depinde de timp, obtinerea functiilor
#  a, H, 𝜺 si calculul curburii k

f = Friedmann
Tmax = 1
tspan = (0.,Tmax)
u0 = [68., 1, 10e-30]
p0 = [1/3., 0.] # w si Λ
prob = ODEProblem(f, u0, tspan, p0)
sol = solve(prob)
plot(sol, vars = (2))
k = sol[2, length(sol)]^2 *( (8*π*G/3)*sol[3, length(sol)]-sol[1, length(sol)]^2 )




f = CSV.read("goodData.csv"; copycols=true)
PolynomFit(t, p) = sqrt.(p[1]*t + p[2]*t.^5)
LineFit(t, p) = p[1]*t .+ p[2]
p0 = [70., 2.]
tdata=f[!, 1]
ydata=float(f[!, 2])

scatter!(tdata, ydata)
fit = curve_fit(LineFit, tdata, ydata, p0)

x = rand(0.1:2.5, 1, 50)
plot(tdata,tdata*fit.param[1] .+fit.param[2], legend=false)
