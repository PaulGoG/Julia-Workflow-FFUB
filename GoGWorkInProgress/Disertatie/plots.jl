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
function Modify_plot(plt::Plots.Plot, xaxisname, yaxisname, xaxislims::Tuple, xaxisscale, yaxislims::Tuple, yaxisscale, plot_title)
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
function Modify_plot(plt::Plots.Plot)
    plot!(plt, minorgrid = true, framestyle = :box)
end
function Plot_textbox(plt::Plots.Plot, x, y, text)
    annotate!(plt, x, y, text)
end
function Plot_legend_attributes(plt::Plots.Plot, lposition)
    plot!(plt, legend_position = lposition)
end
function Process_plot(plt::Plots.Plot, filename::String, fissionant_nucleus_identifier::String)
    savefig(plt, "plots/$(fissionant_nucleus_identifier)_$(filename).png")
    println("*plotting $filename done!")
end
#####
println("*executing data plots")

if secondary_output_Yield == "YES"

end
if secondary_output_ν == "YES"
    plot_surface_ν_A_TKE = Plot_surface(
        DataFrame(x = ν_A_TKE.A, y = ν_A_TKE.TKE, z = ν_A_TKE.Value),
        5, 5, Int(round(10/TKE_step)), 10, (120, 30), 
        "ν(A,TKE)", (minimum(ν_A_TKE.Value), maximum(ν_A_TKE.Value)+0.5), :identity
    )
    display(plot_surface_ν_A_TKE)
    if secondary_output_Ap == "YES"
        plot_surface_y_Ap_Z = Plot_surface(
            DataFrame(x = y_Ap_Z.A, y = y_Ap_Z.Z, z = y_Ap_Z.Value),
            10, 10, 2, 10, (120, 30),
            "Y(Aₚ,Z)", (minimum(y_Ap_Z.Value), maximum(y_Ap_Z.Value)), :identity
        )
        display(plot_surface_y_Ap_Z)
    end
    if secondary_output_Tₖ == "YES"

    end
    if secondary_output_avg_εₖ == "YES"

    end
    if secondary_output_Eᵣ == "YES"

    end
end
if secondary_output_T == "YES"

end
if secondary_output_avg_ε == "YES"

end
if secondary_output_TXE_Q == "YES"
    gr(size = plots_resolution, dpi=300)
    avg_Q = Average_value(Q_AH, y_A, A_H_range)
    plot_Q_AH = Scatter_data(Q_AH.Argument, Q_AH.Value, "", :red, 5, :circle)
    Plot_data(plot_Q_AH, Q_AH.Argument, Q_AH.Value, "", :red)
    Modify_plot(plot_Q_AH)
    Modify_plot(
        plot_Q_AH, L"\mathrm{A_H}", "Q [MeV]", 
        (minimum(Q_AH.Argument), maximum(Q_AH.Argument)), :identity, 
        (minimum(Q_AH.Value)*0.95, maximum(Q_AH.Value)*1.05), :identity, ""
    )
    Plot_textbox(plot_Q_AH, maximum(Q_AH.Argument)*0.95, maximum(Q_AH.Value)*1.025, "<Q> = $(round(avg_Q, digits = 3))")
    xticks!(plot_Q_AH, minimum(Q_AH.Argument):5:maximum(Q_AH.Argument))
    Process_plot(plot_Q_AH, "Q_AH", fissionant_nucleus_identifier)

    avg_TXE = Average_value(txe_AH, y_A, A_H_range)
    plot_TXE_AH = Scatter_data(txe_AH.Argument, txe_AH.Value, "", :red, 5, :circle)
    Plot_data(plot_TXE_AH, txe_AH.Argument, txe_AH.Value, "", :red)
    Modify_plot(plot_TXE_AH)
    Modify_plot(
        plot_TXE_AH, L"\mathrm{A_H}", "TXE [MeV]", 
        (minimum(txe_AH.Argument), maximum(txe_AH.Argument)), :identity, 
        (minimum(txe_AH.Value)*0.95, maximum(txe_AH.Value)*1.05), :identity, ""
    )
    Plot_textbox(plot_TXE_AH, maximum(txe_AH.Argument)*0.85, maximum(txe_AH.Value)*0.8, "<TXE> = $(round(avg_TXE, digits = 3))")
    Process_plot(plot_TXE_AH, "TXE_AH", fissionant_nucleus_identifier)
end
if secondary_output_E_excitation == "YES"

end
if neutron_spectrum == "YES"
    plot_neutron_spectrum = Plot_data(n_E.E, Ratio_to_Maxwellian, "", :red)
    Modify_plot(
        plot_neutron_spectrum, "E [MeV]", "Nuetron spectrum 1/MeV", (first(n_E.E), last(n_E.E)),
        :identity, (minimum(n_E.Value), maximum(n_E.Value)), :log10, "Neutron spectrum, linear scale"
    )
end