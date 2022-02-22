using LsqFit
using Plots
plotlyjs()

g = collect(0.0:2:8.0)
ε = [4880.0, 4080, 3850, 2340, 1850]
σ = [500.0, 700, 970, 1100, 1740]
dε = [ε[1] - ε[i] for i in 2:length(ε)]
x = [g[i] for i in 2:length(g)]

liniarizat(t, p) = p[1] .+ p[2].*t
p0 = [1.0, 1.0]
fit = curve_fit(liniarizat, x, dε, p0)

A = first(fit.param)
B = last(fit.param)
σ_ε = [sqrt(σ[1]^2 + σ[i]^2) for i in 2:length(σ)]

plt = scatter(x, dε, yerr = σ_ε,
    title = "Graficul atenuării liniare a particulei α prin mylar",
    titlefont = (12, "times"),
    xlabel = "x (μm)",
    ylabel = "dε (KeV)",
    label = "mylar",
    framestyle = :box,
    color = :black)
plt = plot!(x, y -> A + B * y, label = false)
display(plt)