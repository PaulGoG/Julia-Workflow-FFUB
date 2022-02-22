using LsqFit
using Plots
plotlyjs()

t_a = collect(100:100.0:500)
N_Mg = [1459.0, 1018, 1175, 953, 829]
N_Al = [2985.0, 1728, 1116, 632, 402]
N_Mg = log.(N_Mg[1]./N_Mg)
N_Al = log.(N_Al[1]./N_Al)

liniarizat(t, p) = p[1] .+ p[2].*t
p0 = [1.0, 1.0]
fit_Mg = curve_fit(liniarizat, t_a, N_Mg, p0)
fit_Al = curve_fit(liniarizat, t_a, N_Al, p0)

lambda_Mg = last(fit_Mg.param)
lambda_Al = last(fit_Al.param)
A_Mg = first(fit_Mg.param)
A_Al = first(fit_Al.param)
T_Mg = (log(2)/lambda_Mg)/60
T_Al = (log(2)/lambda_Al)/60

plt = scatter(t_a, N_Mg,
    title = "Legea Schweidler liniarizată",
    titlefont = (12, "times"),
    xlabel = "t_achizitie (s)",
    ylabel = "ln(N₀/N)",
    label = "Mg",
    framestyle = :box,
    color = :black)
plt = scatter!(t_a, N_Al, label = "Al", color = :red)
plt = plot!(t_a, y -> A_Mg + lambda_Mg * y, label = false)
plt = plot!(t_a, y -> A_Al + lambda_Al * y, label = false)
println(T_Mg)
println(T_Al)
display(plt)
#savefig("3_grafic")