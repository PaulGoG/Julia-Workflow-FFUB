using LsqFit
using Plots
plotlyjs()

λ = [7091.99, 7081.88, 6907.16, 6716.17, 5790.65, 5789.66, 5460.74, 4916.04, 4358.35, 4077.81, 4046.56]
x = [2.0, 7, 15, 17, 35, 36, 60, 85, 120, 180, 200]

Cauchy(t, p) = p[1] .+ p[2]./t.^2
p0 = [1.0, 1.0]
fit = curve_fit(Cauchy, λ, x, p0)

A = first(fit.param)
B = last(fit.param)

plt = scatter(λ, x,
    title = "Etalonarea spectroscopului cu lampa de Hg",
    titlefont = (12, "times"),
    xlabel = "λ (Å)",
    ylabel = "x (div)",
    label = false,
    framestyle = :box,
    color = :black
)

x_plot = collect(minimum(λ)-300:100.0:maximum(λ)+300)
plt = plot!(x_plot, y -> A + B / y^2, label = false)
#savefig(plt, "Calibrare.png")