using QuadGK, Plots, DataFrames, LaTeXStrings

cd(@__DIR__)

Rₐ = 50e-3    #Raza amplasament in km
Pᵢ = 1e-9     #Probabilitatea de pierdere a controlului per km de zbor
Nc = 7e4      #Nr. de zboruri anuale pe ruta analizata
g = 0.23      #Densitatea de probabilitate unitara de prabusire a avionului dupa pierderea controlului in (r, r-dr)

y₀_Min = 5       #Distanta minima dintre centrul amplasamentului si ruta de zbor in km
y₀_Max = 100      #Distanta maxima dintre centrul amplasamentului si ruta de zbor in km
x₀ = 200          #Lungimea rutei de zbor pe care se realizeaza integrarea in km

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
        push!(freqAnualaPrabusire.fR, Pᵢ *Nc *π *R^2 *integralaRadiala)
        push!(freqAnualaPrabusire.fRU, Pᵢ *Nc *π *R^2 *integralaRadialaUnghiulara)
    end
end


#Reprezentari grafice
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
        c = color_scale,
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
function Bar_data(x, y, plot_label, color_bars, Δx)
    plt = bar(
        x, y, label = plot_label, 
        linetype = :steppre, fillalpha = 0, linecolor = color_bars, bar_width = Δx
    )
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
    println("*plotting $filename done!")
end

gr(size = (900, 900), dpi=600)
xData = freqAnualaPrabusire.y₀
yDataR = freqAnualaPrabusire.fR
yDataRU = freqAnualaPrabusire.fRU

plotFreq2DLin = Plot_data(xData, yDataR, "FDP radiala", :red)
Plot_data(plotFreq2DLin, xData, yDataRU, "FDP radiala si unghiulara", :blue)
Modify_plot(plotFreq2DLin)
Modify_plot(
    plotFreq2DLin, "y₀", "Frecventa anuala de prabusire", 
    (minimum(xData), maximum(xData)), :identity, 
    (0.0, maximum(yDataR)*1.1), :identity, ""
)

Process_plot(plot_y_A_lin, "y_A_lin", fissionant_nucleus_identifier)

####

plot_nSpectrum_ratio_Maxwellian = Plot_data(n_E_SL.E, Ratio_to_Maxwellian_SL, "", :red)
Modify_plot(
    plot_nSpectrum_ratio_Maxwellian, "E [MeV]", "Ratio to Maxwellian spectrum", (first(n_E_SL.E), 10.0),
    :log10, (0.5, 1.5), :identity, "SL Neutron spectrum as ratio to Maxwellian"
)
hline!(plot_nSpectrum_ratio_Maxwellian, [1.0], linestyle = :dashdot, color = :black, label = latexstring("\$\\mathrm{T_M}\$ = $(round(T_M_eq_SL, digits = 2)) MeV"))
Modify_plot(plot_nSpectrum_ratio_Maxwellian)
xticks!(plot_nSpectrum_ratio_Maxwellian, [10.0^i for i in -10:10])
Process_plot(plot_nSpectrum_ratio_Maxwellian, "nSpectrum_ratio_Maxwellian", fissionant_nucleus_identifier)

####

plotlyjs(size = (16, 9) .* 90, dpi=600)
plot_surface_ν_A_TKE = Plot_surface(
    DataFrame(x = ν_A_TKE.A, y = ν_A_TKE.TKE, z = ν_A_TKE.Value),
    10, 10, 20, Int(round(10/TKE_step)), (120, 0), 
    "ν(A,TKE)", (0.0, maximum(ν_A_TKE.Value) + 0.5), :identity,
    :turbo
)
Modify_plot(plot_surface_ν_A_TKE, "A", "TKE", "")
Modify_plot(plot_surface_ν_A_TKE)
display(plot_surface_ν_A_TKE)

####