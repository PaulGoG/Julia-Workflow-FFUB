include("preamble.jl")
include("readData.jl")
include("buildTimeseries.jl")
include("parameterOptimization.jl")
include("optimizeDelta.jl")
include("plotFigs.jl")

if δ_optimization
    δ = optimumDelta(datapoints, δ)
end

println("*parameter optimization for δ = $δ inbound...\n")
timeseries = wtsProbabilities(waiting_time_series(datapoints, δ))
α, xₘ, d_KS = powerlawOptimize(timeseries)

println("α = $(round(α, digits = noDigits))")
println("xₘ = $(round(xₘ, digits = noDigits))")
println("KS = $(round(d_KS, digits = noDigits))\n")
println("*executing plot figures...")

gr(size = (1000, 1000), dpi = 600)
P_PDF(x, α, xₘ) = (α - 1)/xₘ * (xₘ/x)^α

#Main 2D plot of PDF together with powerlaw fit of data at fixed (given or optimal) δ value
xFScale = @view timeseries.time[timeseries.time .>= xₘ]
yFScale = P_PDF.(xFScale, α, xₘ)
fitLabel =  string("KS = $(round(d_KS, digits = noDigits)); ",
            "α = $(round(α, digits = noDigits)); ",
            "xₘ = $(round(xₘ, digits = noDigits)); ",
            "δ = $δ")

plotPDF_2D = Scatter_data(timeseries.time, timeseries.pdf, "", :blue, 3, :xcross)
Plot_data(plotPDF_2D, xFScale, yFScale, fitLabel, :red)
Modify_plot(plotPDF_2D)
Modify_plot(
    plotPDF_2D, "time ($datafileTime)", "Pₖ", 
    (0.1, maximum(timeseries.time)*1e1), :log10, 
    (minimum(timeseries.pdf)*1e-1, 1.0), :log10, "Probability Density Function of $datafileQuantity"
)
Plot_textbox(plotPDF_2D, 2.5, minimum(timeseries.pdf)*2e-1, 
            "Data source: $datafileName")
xticks!(plotPDF_2D, [10^i for i in 0:10])
yticks!(plotPDF_2D, [10.0^i for i in -10:0])
Process_plot(plotPDF_2D, "PDF_$(datafileName)_col$(datafileCol)")
println("*done plotting PDF_$(datafileName)_col$(datafileCol)\n")

#2D plot of KS(δ) at optimal (α, xₘ) values for each δ
δ_range = δ*1e-2:δ*1e-1:δ*2
KS_vals = Float64[]
for δ_vals in δ_range
    println("plotting KS for δ = $δ_vals")
    ts = wtsProbabilities(waiting_time_series(datapoints, δ_vals))
    push!(KS_vals, powerlawOptimize(ts)[3])
end

plotKSδ_2D = Scatter_data(δ_range, KS_vals, "", :red, 5, :circle)
Scatter_data(plotKSδ_2D, [δ], [d_KS], "δ = $δ", :black, 8, :star)
Modify_plot(plotKSδ_2D)
Modify_plot(
    plotKSδ_2D, "δ", "KS", 
    (minimum(δ_range) - δ*1e-1, maximum(δ_range) + δ*1e-1), :identity, 
    (minimum(KS_vals)*0.9, maximum(KS_vals)*1.1), :identity, "Data source: $datafileName"
)
display(plotKSδ_2D)
Process_plot(plotKSδ_2D, "KSdelta_$(datafileName)_col$(datafileCol)")
println("\n*done plotting KSdelta_$(datafileName)_col$(datafileCol)")

#3D plots (surface and heatmap) of KS values at various (α, xₘ) for fixed (given or optimal) δ value
if triDimPlots
    plotlyjs(size = (16, 9) .* 90, dpi = 600)

    KSDF = DataFrame(α = Float64[], xₘ = Float64[], KS = Float64[])
    for α in 0.9*α:α*1e-2:α*1.1
        for xₘ in 0.0:xₘ*0.1:xₘ*10
            push!(KSDF.α, α)
            push!(KSDF.xₘ, xₘ)
            push!(KSDF.KS, lossKolmogorovSmirnov((α, xₘ), timeseries.time, timeseries.cdf))
        end
    end

    plotKS_surface = Plot_surface(
        DataFrame(x = KSDF.α, y = KSDF.xₘ, z = KSDF.KS),
        1, α*(0.9 + 1e-2), 10, xₘ, (120, 0), 
        "KS distance", (minimum(KSDF.KS), maximum(KSDF.KS)), :identity,
        :turbo
    )
    Modify_plot(plotKS_surface, "α", "xₘ", "")
    Modify_plot(plotKS_surface)
    display(plotKS_surface)

    plotKS_heatmap = Plot_heatmap(
        DataFrame(x = KSDF.α, y = KSDF.xₘ, z = KSDF.KS),
        1, α*(0.9 + 1e-2), 10, xₘ,
        "KS distance", :turbo
    )
    Modify_plot(plotKS_heatmap, "α", "xₘ", "")
    Modify_plot(plotKS_heatmap)
    display(plotKS_heatmap)
end


println("*ending program execution at $(Dates.format(now(), "HH:MM:SS"))")
println("*program execution succesful!")