#Function bodies and function calls for plotting data
#####
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
function Plot_surface(q_x_y::DataFrame, tick_size_xaxis::Int, tick_roundness_xaxis, tick_size_yaxis::Int, tick_roundness_yaxis, camera_angle::Tuple, zaxisname, zaxislims::Tuple, zaxisscale)
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
        zlabel = zaxisname,
        zlims = zaxislims,
        zscale = zaxisscale
    )
    return plt
end
function Modify_plot(plt::Plots.Plot, xaxisname, yaxisname, xaxislims::Tuple, xaxisscale, yaxislims::Tuple, yaxisscale, plot_title, DPI::Int)
    plot!(
        plt,
        xlabel = xaxisname,
        ylabel = yaxisname,
        xlims = xaxislims,
        ylims = yaxislims,
        xscale = xaxisscale,
        yscale = yaxisscale,
        title = plot_title,
        dpi = DPI,
        minorgrid = true,
        framestyle = :box
    )
end
function Plot_data(x, y, plot_label, plot_color)
    plt = plot(x, y, label = plot_label, color = plot_color)
    return plt
end
function Plot_data(x, y, σ, plot_label, plot_color)
    plt = plot(x, y, ribbon = σ, label = plot_label, color = plot_color)
    return plt
end
function Scatter_data(x, y, plot_label, marker_color, marker_size, marker_shape)
    plt = scatter(
        x, y, label = plot_label, 
        markercolor = marker_color, markersize = marker_size, markershape = marker_shape
    )
    return plt
end
function Scatter_data(x, y, σ, plot_label, marker_color, marker_size, marker_shape)
    plt = scatter(
        x, y, yerror = σ, label = plot_label, 
        markercolor = marker_color, markersize = marker_size, markershape = marker_shape
    )
    return plt
end
function Plot_data(plt::Plots.Plot, x, y, plot_label, plot_color)
    plot!(plt, x, y, label = plot_label, color = plot_color)
end
function Plot_data(plt::Plots.Plot, x, y, σ, plot_label, plot_color)
    plot!(plt, x, y, ribbon = σ, label = plot_label, color = plot_color)
end
function Scatter_data(plt::Plots.Plot, x, y, plot_label, marker_color, marker_size, marker_shape)
    scatter!(
        plt,
        x, y, label = plot_label, 
        markercolor = marker_color, markersize = marker_size, markershape = marker_shape
    )
end
function Scatter_data(plt::Plots.Plot, x, y, σ, plot_label, marker_color, marker_size, marker_shape)
    scatter!(
        plt,
        x, y, yerror = σ, label = plot_label, 
        markercolor = marker_color, markersize = marker_size, markershape = marker_shape
    )
end
function Plot_textbox(plt::Plots.Plot, x, y, text)
    annotate!(plt, x, y, text)
end
function Plot_legend_attributes(plt::Plots.Plot, lposition)
    plot!(legend_position = lposition)
end
function Display_plot(plt::Plots.Plot, filename::String)
    display(plt)
    #savefig("plots/$filename.png")
end
#####
#=
plot_surface_ν_A_TKE = Plot_surface(
    DataFrame(x = ν_A_TKE.A, y = ν_A_TKE.TKE, z = ν_A_TKE.Value),
    10, 10, Int(round(10/TKE_step)), 10, (120, 30), 
    "ν", (minimum(ν_A_TKE.Value), maximum(ν_A_TKE.Value)), :identity
    )
display(plot_surface_ν_A_TKE)

plot_surface_y_Ap_Z = Plot_surface(
    DataFrame(x = y_Ap_Z.A, y = y_Ap_Z.Z, z = y_Ap_Z.Value),
    10, 10, 2, 10, (120, 30),
    "Y(Aₚ,Z)", (minimum(y_Ap_Z.Value), maximum(y_Ap_Z.Value)), :identity
    )
display(plot_surface_y_Ap_Z)
=#

plot_neutron_spectrum = Plot_data(n_E.E, n_E.Value, "", :red)
plot_neutron_spectrum = Modify_plot(
    plot_neutron_spectrum, "E [MeV]", "N(E)", (first(n_E.E), last(n_E.E)),
    :identity, (first(n_E.Value), last(n_E.Value)), :log10, "Neutron spectrum in SL", 300 
)
Display_plot(plot_neutron_spectrum, "randomfilename")