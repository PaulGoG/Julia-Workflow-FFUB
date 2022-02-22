using LsqFit
using Plots
plotlyjs()

Energii = [511.0, 1173.0, 1274.0, 1332.0] #KeV
Canale = [77.0, 172.0, 185.0, 195.0]

liniarizat(t, p) = p[1] .+ p[2].*t
p0 = [1.0, 1.0]
# Energia = A + B*Canal, calibrarea energetica a detectorului
fit_calibrare = curve_fit(liniarizat, Canale, Energii, p0)
A = first(fit_calibrare.param)
B = last(fit_calibrare.param)
Energie_Pic_Studiat = round(A + 67*B, digits = 2);

plt = scatter(Canale, Energii,
    title = "Calibrarea energetică a detectorului",
    titlefont = (12, "times"),
    xlabel = "Nr. canal",
    ylabel = "ε (KeV)",
    legend = false,
    framestyle = :box,
    color = :black)
plt = plot!(Canale, y -> A + B*y, color = :magenta)
display(plt)
#savefig("7_grafic1")

Timpi = collect(520.0:260:1040)
N_exp = [3691.0, 3275.0, 2827.0, 2487.0]
N = [log.(N_exp[1]./N_exp[i]) for i in 2:length(N_exp)]
M = [N_exp[i] for i in 2:length(N_exp)]
wt = M
fit_fotopic = curve_fit(liniarizat, Timpi, N, wt, p0)
lambda = last(fit_fotopic.param)
T_inj = (log(2)/lambda)/60

plt = scatter(Timpi, N, yerr = 1 ./sqrt.(M),
    title = "Legea Schweidler liniarizată",
    titlefont = (12, "times"),
    xlabel = "t_achizitie (s)",
    ylabel = "ln(N₀/N)",
    legend = false,
    framestyle = :box,
    color = :purple)
plt = plot!(Timpi, y -> first(fit_fotopic.param) + lambda * y, color = :red)
display(plt)
#savefig("7_grafic2")

println("Energie pic studiat = ", Energie_Pic_Studiat, " KeV")
σ_λ = last(stderror(fit_fotopic));
println("λ = ", lambda," +- ",σ_λ)
println("T1/2 = ", T_inj," +- ", ((log(2)/lambda^2)/60) *σ_λ)