using QuadGK, Plots, DataFrames, LaTeXStrings

cd(@__DIR__)

Rₐ = 50       #Raza amplasament in m
Pᵢ = 1e-9     #Probabilitatea de pierdere a controlului per km de zbor
Nc = 7e4      #Nr. de zboruri anuale pe ruta analizata
g = 0.23      #Densitatea de probabilitate unitara de prabusire a avionului dupa pierderea controlului in (r, r-dr)

y₀_Min = 5       #Distanta minima dintre centrul amplasamentului si ruta de zbor in km
y₀_Max = 50      #Distanta maxima dintre centrul amplasamentului si ruta de zbor in km
x₀ = 200         #Lungimea rutei de zbor pe care se realizeaza integrarea in km

freqAnualaPrabusire = DataFrame(R = Float64[], y₀ = Float64[], fR = Float64[], fRU = Float64[])

for y₀ in y₀_Min:y₀_Max
    function integrandRadial(x)
        if x > 0
            return x/sqrt(x^2 + y₀^2) *exp(-g * sqrt(x^2 + y₀^2))
        else
            return 0.0
        end
    end
    function integrandRadialUnghiular(x)
        if x > 0
            return x*y₀/(x^2 + y₀^2) *exp(-g * sqrt(x^2 + y₀^2))
        else
            return 0.0
        end
    end
    integralaRadiala = first(quadgk(integrandRadial, 0.0, x₀))
    integralaRadialaUnghiulara = g/2 *first(quadgk(integrandRadialUnghiular, 0.0, x₀))
    for R in Rₐ*1e-1:Rₐ*1e-1:Rₐ*10
        push!(freqAnualaPrabusire.R, R)
        push!(freqAnualaPrabusire.y₀, y₀)
        push!(freqAnualaPrabusire.fR, Pᵢ *Nc *π *(R*1e-3)^2 *integralaRadiala)
        push!(freqAnualaPrabusire.fRU, Pᵢ *Nc *π *(R*1e-3)^2 *integralaRadialaUnghiulara)
    end
end
################################
################################
################################
function Grid_builder(q_x_y::DataFrame)
    x = sort(unique(q_x_y.x))
    y = sort(unique(q_x_y.y))
    grid = zeros(length(x), length(y))
    for index_x in eachindex(x)
        for index_y in eachindex(y)
            if isassigned(q_x_y.z[(q_x_y.x .== x[index_x]) .& (q_x_y.y .== y[index_y])], 1)
                grid[index_x, index_y] = q_x_y.z[(q_x_y.x .== x[index_x]) .& (q_x_y.y .== y[index_y])][1] 
            end
        end
    end
    return grid
end
function Plot_surface(q_x_y::DataFrame, tick_size_xaxis::Int, tick_roundness_xaxis, tick_size_yaxis::Int, tick_roundness_yaxis, camera_angle::Tuple, zaxisname, zaxislims::Tuple, zaxisscale, color_scale)
    grid_q_x_y = Grid_builder(q_x_y)

    x = sort(unique(q_x_y.x))
    first_index_x = findfirst(x -> round(x/tick_roundness_xaxis) == x/tick_roundness_xaxis, x)
    last_index_x = findlast(x -> round(x/tick_roundness_xaxis) == x/tick_roundness_xaxis, x)
    index_x_range = first_index_x:tick_size_xaxis:last_index_x

    y = sort(unique(q_x_y.y))
    first_index_y = findfirst(x -> round(x/tick_roundness_yaxis) == x/tick_roundness_yaxis, y)
    last_index_y = findlast(x -> round(x/tick_roundness_yaxis) == x/tick_roundness_yaxis, y)
    index_y_range = first_index_y:tick_size_yaxis:last_index_y

    plt = surface(
        grid_q_x_y', 
        xticks = (collect(index_x_range), string.(Int.(round.(x[index_x_range])))), 
        yticks = (collect(index_y_range), string.(Int.(round.(y[index_y_range])))),
        camera = camera_angle,
        c = cgrad(color_scale, scale = zaxisscale),
        zlabel = zaxisname,
        zlims = zaxislims,
        zscale = zaxisscale,
        colorbar_title = zaxisname
    )
    return plt
end
function Plot_heatmap(q_x_y::DataFrame, tick_size_xaxis::Int, tick_roundness_xaxis, 
    tick_size_yaxis::Int, tick_roundness_yaxis, zaxisname, color_scale)
    grid_q_x_y = Grid_builder(q_x_y)

    x = sort(unique(q_x_y.x))
    first_index_x = findfirst(x -> round(x/tick_roundness_xaxis) == x/tick_roundness_xaxis, x)
    last_index_x = findlast(x -> round(x/tick_roundness_xaxis) == x/tick_roundness_xaxis, x)
    index_x_range = first_index_x:tick_size_xaxis:last_index_x

    y = sort(unique(q_x_y.y))
    first_index_y = findfirst(x -> round(x/tick_roundness_yaxis) == x/tick_roundness_yaxis, y)
    last_index_y = findlast(x -> round(x/tick_roundness_yaxis) == x/tick_roundness_yaxis, y)
    index_y_range = first_index_y:tick_size_yaxis:last_index_y

    plt = heatmap(
        grid_q_x_y', 
        xticks = (collect(index_x_range), string.(Int.(round.(x[index_x_range])))), 
        yticks = (collect(index_y_range), string.(Int.(round.(y[index_y_range])))),
        c = color_scale,
        colorbar_title = zaxisname,
    )
    return plt
end
function Plot_data(x, y, plot_label, plot_color)
    plt = plot(x, y, label = plot_label, color = plot_color)
    return plt
end
function Plot_data(plt::Plots.Plot, x, y, plot_label, plot_color)
    plot!(plt, x, y, label = plot_label, color = plot_color)
end
function Modify_plot(plt::Plots.Plot, xaxisname, yaxisname, xaxislims::Tuple, 
    xaxisscale, yaxislims::Tuple, yaxisscale, plot_title)
    plot!(
        plt,
        xlabel = xaxisname,
        ylabel = yaxisname,
        xlims = xaxislims,
        ylims = yaxislims,
        xscale = xaxisscale,
        yscale = yaxisscale,
        title = plot_title,
    )
end
function Modify_plot(plt::Plots.Plot, xaxisname, yaxisname, plot_title)
    plot!(
        plt,
        xlabel = xaxisname,
        ylabel = yaxisname,
        title = plot_title,
    )
end
function Modify_plot(plt::Plots.Plot)
    plot!(plt, minorgrid = true, framestyle = :box)
end
function Plot_textbox(plt::Plots.Plot, x, y, text)
    annotate!(plt, x, y, text)
end
function Plot_legend_attributes(plt::Plots.Plot, lposition)
    plot!(plt, legend_position = lposition)
end
function Process_plot(plt::Plots.Plot, filename::String)
    savefig(plt, "$(filename).png")
end

gr(size = (900, 900), dpi=600)
xData = freqAnualaPrabusire.y₀[freqAnualaPrabusire.R .== Rₐ]
yDataR = freqAnualaPrabusire.fR[freqAnualaPrabusire.R .== Rₐ]
yDataRU = freqAnualaPrabusire.fRU[freqAnualaPrabusire.R .== Rₐ]

plotFreq2DLin = Plot_data(xData, yDataR, L"\mathrm{f_R}(r)", :red)
Plot_data(plotFreq2DLin, xData, yDataRU, L"\mathrm{f_R}(r, \theta)", :blue)
Modify_plot(plotFreq2DLin)
Modify_plot(
    plotFreq2DLin, "y₀ (km)", "Frecventa anuala de prabusire", 
    (minimum(xData), maximum(xData)), :identity, 
    (0.0, maximum(yDataR)), :identity, ""
)
xticks!(plotFreq2DLin, y₀_Min:5:y₀_Max)
yticks!(plotFreq2DLin, [i*1e-8 for i in 0:10:100])
display(plotFreq2DLin)
#Process_plot(plotFreq2DLin, "FrecventaPrabusire2Dliniar")

plotFreq2DLog = Plot_data(xData, yDataR, L"\mathrm{f_R}(r)", :red)
Plot_data(plotFreq2DLog, xData, yDataRU, L"\mathrm{f_R}(r, \theta)", :blue)
Modify_plot(plotFreq2DLog)
Modify_plot(
    plotFreq2DLog, "y₀ (km)", "Frecventa anuala de prabusire", 
    (minimum(xData), maximum(xData)), :identity, 
    (minimum(yDataRU)*1e-1, 1.0), :log10, ""
)
hline!(plotFreq2DLog, [1e-5], linestyle = :dashdot, color = :black, label = "Accident baza de proiectare de tip A")
hline!(plotFreq2DLog, [1e-7], linestyle = :dash, color = :black, label = "Accident baza de proiectare de tip B")
yticks!(plotFreq2DLog, [10.0^i for i in -20:0])
xticks!(plotFreq2DLog, y₀_Min:5:y₀_Max)
display(plotFreq2DLog)
#Process_plot(plotFreq2DLog, "FrecventaPrabusire2Dlogaritmic")

plotlyjs(size = (16, 9) .* 90, dpi=600)

plotFreq3DSuprafata = Plot_surface(
    DataFrame(x = freqAnualaPrabusire.R, y = freqAnualaPrabusire.y₀, z = freqAnualaPrabusire.fRU),
    10, 50, 5, 5, (120, 0), 
    "Freq", (minimum(freqAnualaPrabusire.fRU), maximum(freqAnualaPrabusire.fRU)), :identity,
    :turbo
)
Modify_plot(plotFreq3DSuprafata, "R (m)", "y₀ (km)", "")
Modify_plot(plotFreq3DSuprafata)
display(plotFreq3DSuprafata)

plotFreq3DHeatmap = Plot_heatmap(
    DataFrame(x = freqAnualaPrabusire.R, y = freqAnualaPrabusire.y₀, z = freqAnualaPrabusire.fRU),
    10, 50, 5, 5, 
    "Freq", :turbo
)
Modify_plot(plotFreq3DHeatmap, "R (m)", "y₀ (km)", "")
Modify_plot(plotFreq3DHeatmap)
display(plotFreq3DHeatmap)