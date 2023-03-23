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
function Plot_heatmap(q_x_y::DataFrame, tick_size_xaxis::Int, tick_roundness_xaxis, tick_size_yaxis::Int, tick_roundness_yaxis, zaxisname, color_scale)
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
function Bar_data(x, y, plot_label, color_bars, Δx)
    plt = bar(
        x, y, label = plot_label, 
        linetype = :steppre, fillalpha = 0, linecolor = color_bars, bar_width = Δx
    )
    return plt
end
function Bar_data(x, y, σ, plot_label, color_bars, Δx)
    plt = bar(
        x, y, yerror = σ, label = plot_label, 
        linetype = :steppre, fillalpha = 0, linecolor = color_bars, bar_width = Δx
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
function Bar_data(plt::Plots.Plot, x, y, plot_label, color_bars, Δx)
    bar!(
        plt,
        x, y, label = plot_label, 
        linetype = :steppre, fillalpha = 0, linecolor = color_bars, bar_width = Δx
    )
end
function Bar_data(plt::Plots.Plot, x, y, σ, plot_label, color_bars, Δx)
    bar!(
        plt,
        x, y, yerror = σ, label = plot_label, 
        linetype = :steppre, fillalpha = 0, linecolor = color_bars, bar_width = Δx
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
function Process_plot(plt::Plots.Plot, filename::String, fissionant_nucleus_identifier::String)
    savefig(plt, "plots/$(fissionant_nucleus_identifier)_$(filename).png")
    println("*plotting $filename done!")
end
#####
println("*executing data plots")

if secondary_output_Yield == "YES"
    gr(size = plots_resolution, dpi=300)

    avg_A_L = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .<= A_H_min])
    avg_A_H = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .>= A_H_min])
    plot_y_A = Scatter_data(y_A.Argument, y_A.Value, "", :red, 5, :circle)
    Plot_data(plot_y_A, y_A.Argument, y_A.Value, "", :red)
    Modify_plot(plot_y_A)
    Modify_plot(
        plot_y_A, "A", "Yield %", 
        (minimum(y_A.Argument), maximum(y_A.Argument)), :identity, 
        (minimum(y_A.Value)*0.95, maximum(y_A.Value)*1.05), :identity, ""
    )
    Plot_textbox(plot_y_A, maximum(y_A.Argument)*0.95, maximum(y_A.Value)*1.025, latexstring("\$\\mathrm{<A_L>}\$ = $(round(avg_A_L[1], digits = 3))"))
    Plot_textbox(plot_y_A, maximum(y_A.Argument)*0.95, maximum(y_A.Value)*1.015, latexstring("\$\\mathrm{<A_H>}\$ = $(round(avg_A_H[1], digits = 3))"))
    xticks!(plot_y_A, 10 *div(minimum(y_A.Argument), 10):5:maximum(y_A.Argument))
    Process_plot(plot_y_A, "y_A_lin", fissionant_nucleus_identifier)
    Modify_plot(
        plot_y_A, "A", "Yield %", 
        (minimum(y_A.Argument), maximum(y_A.Argument)), :identity, 
        (minimum(y_A.Value)*0.95, maximum(y_A.Value)*1.05), :log10, ""
    )
    Process_plot(plot_y_A, "y_A_log", fissionant_nucleus_identifier)

    
end
if secondary_output_ν == "YES"
    plotlyjs(size = (16, 9) .* 90, dpi=300)

    plot_surface_ν_A_TKE = Plot_surface(
        DataFrame(x = ν_A_TKE.A, y = ν_A_TKE.TKE, z = ν_A_TKE.Value),
        10, 5, 20, Int(round(10/TKE_step)), (120, 0), 
        "ν(A,TKE)", (minimum(ν_A_TKE.Value), maximum(ν_A_TKE.Value)+0.5), :identity,
        :turbo
    )
    Modify_plot(plot_surface_ν_A_TKE, "A", "TKE", "")
    Modify_plot(plot_surface_ν_A_TKE)
    #display(plot_surface_ν_A_TKE)

    gr(size = plots_resolution, dpi=300)

    avg_ν_L = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_L_range)
    avg_ν_H = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_H_range)
    plot_ν_A = Scatter_data(ν_A.Argument, ν_A.Value, "", :red, 5, :circle)
    Plot_data(plot_ν_A, ν_A.Argument, ν_A.Value, "", :red)
    Modify_plot(plot_ν_A)
    Modify_plot(
        plot_ν_A, "A", "ν", 
        (minimum(ν_A.Argument), maximum(ν_A.Argument)), :identity, 
        (minimum(ν_A.Value)*0.95, maximum(ν_A.Value)*1.05), :identity, ""
    )
    Plot_textbox(plot_ν_A, maximum(ν_A.Argument)*0.95, maximum(ν_A.Value)*1.025, latexstring("\$\\mathrm{<\\nu_L>}\$ = $(round(avg_ν_L, digits = 3))"))
    Plot_textbox(plot_ν_A, maximum(ν_A.Argument)*0.95, maximum(ν_A.Value)*1.015, latexstring("\$\\mathrm{<\\nu_H>}\$ = $(round(avg_ν_H, digits = 3))"))
    xticks!(plot_ν_A, 10 *div(minimum(ν_A.Argument), 10):5:maximum(ν_A.Argument))
    Process_plot(plot_ν_A, "nu_A", fissionant_nucleus_identifier)

    plot_P_ν = Scatter_data(probability_ν.Argument, probability_ν.Value, "", :red, 5, :circle)
    Plot_data(plot_P_ν, probability_ν.Argument, probability_ν.Value, "", :red)
    Modify_plot(plot_P_ν)
    Modify_plot(
        plot_P_ν, "ν", "Probability %", 
        (minimum(probability_ν.Argument), maximum(probability_ν.Argument)), :identity, 
        (minimum(probability_ν.Value)*0.95, maximum(probability_ν.Value)*1.05), :identity, ""
    )
    Process_plot(plot_P_ν, "P_nu", fissionant_nucleus_identifier)

    if secondary_output_Ap == "YES"
        plotlyjs(size = (16, 9) .* 90, dpi=300)

        plot_heatmap_y_Ap_Z = Plot_heatmap(
            DataFrame(x = y_Ap_Z.A, y = y_Ap_Z.Z, z = y_Ap_Z.Value),
            10, 10, 2, 2,
            "Y(Aₚ,Z)", :turbo
        )
        Modify_plot(plot_heatmap_y_Ap_Z, "Aₚ", "Z", "")
        Modify_plot(plot_heatmap_y_Ap_Z)
        #display(plot_heatmap_y_Ap_Z)

        gr(size = plots_resolution, dpi=300)

        avg_Ap_L = Average_yield_argument(y_Ap, y_Ap.Argument[y_Ap.Argument .< Ap_H_min])
        avg_Ap_H = Average_yield_argument(y_Ap, y_Ap.Argument[y_Ap.Argument .>= Ap_H_min])
        plot_y_Ap = Scatter_data(y_Ap.Argument, y_Ap.Value, "", :red, 5, :circle)
        Plot_data(plot_y_Ap, y_Ap.Argument, y_Ap.Value, "", :red)
        Modify_plot(plot_y_Ap)
        Modify_plot(
            plot_y_Ap, "Aₚ", "Yield %", 
            (minimum(y_Ap.Argument), maximum(y_Ap.Argument)), :identity, 
            (minimum(y_Ap.Value)*0.95, maximum(y_Ap.Value)*1.05), :identity, ""
        )
        Plot_textbox(plot_y_Ap, maximum(y_Ap.Argument)*0.95, maximum(y_Ap.Value)*1.025, latexstring("\$\\mathrm{<Ap_L>}\$ = $(round(avg_Ap_L[1], digits = 3))"))
        Plot_textbox(plot_y_Ap, maximum(y_Ap.Argument)*0.95, maximum(y_Ap.Value)*1.015, latexstring("\$\\mathrm{<Ap_H>}\$ = $(round(avg_Ap_H[1], digits = 3))"))
        xticks!(plot_y_Ap, 10 *div(minimum(y_Ap.Argument), 10):5:maximum(y_Ap.Argument))
        Process_plot(plot_y_Ap, "y_Ap_lin", fissionant_nucleus_identifier)
    end
    if secondary_output_Tₖ == "YES"

    end
    if secondary_output_avg_εₖ == "YES"

    end
    if secondary_output_Eᵣ == "YES"

    end
end
if secondary_output_T == "YES"
    gr(size = plots_resolution, dpi=300)

    plot_P_T = Bar_data(probability_T.Argument, probability_T.Value, "", :red, ΔT)
    Modify_plot(plot_P_T)
    Modify_plot(
        plot_P_T, "T [MeV]", "Probability %", 
        (0.0, maximum(probability_T.Argument)), :identity, 
        (0.0, maximum(probability_T.Value)*1.05), :identity, ""
    )
    Process_plot(plot_P_T, "P_T", fissionant_nucleus_identifier)
end
if secondary_output_avg_ε == "YES"
    gr(size = plots_resolution, dpi=300)

    plot_avgE_SCM = Bar_data(probability_avg_ε.Argument, probability_avg_ε.Value, "", :red, Δavg_ε)
    Modify_plot(plot_avgE_SCM)
    Modify_plot(
        plot_avgE_SCM, "<ε> [MeV]", "Probability %", 
        (0.0, maximum(probability_avg_ε.Argument)), :identity, 
        (0.0, maximum(probability_avg_ε.Value)*1.05), :identity, ""
    )
    Process_plot(plot_avgE_SCM, "P_avgE_SCM", fissionant_nucleus_identifier)
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
    xticks!(plot_Q_AH, 10 *div(minimum(Q_AH.Argument), 10):5:maximum(Q_AH.Argument))
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
    Plot_textbox(plot_TXE_AH, maximum(txe_AH.Argument)*0.95, maximum(txe_AH.Value)*1.025, "<TXE> = $(round(avg_TXE, digits = 3))")
    xticks!(plot_TXE_AH, 10 *div(minimum(txe_AH.Argument), 10):5:maximum(txe_AH.Argument))
    Process_plot(plot_TXE_AH, "TXE_AH", fissionant_nucleus_identifier)
end
if secondary_output_E_excitation == "YES"
    gr(size = plots_resolution, dpi=300)

    avg_E_exi = Average_value(E_excitation_A, y_A, A_range)
    plot_P_E = Bar_data(probability_E_excitation.Argument, probability_E_excitation.Value, "", :red, 1.0)
    Modify_plot(plot_P_E)
    Modify_plot(
        plot_P_E, "E* [MeV]", "Probability %", 
        (0.0, 40.0), :identity, 
        (0.0, maximum(probability_E_excitation.Value)*1.05), :identity, ""
    )
    Plot_textbox(plot_P_E, 40.0*0.85, maximum(probability_E_excitation.Value)*1.02, "<E*> = $(round(avg_E_exi, digits = 3))")
    xticks!(plot_P_E, 10 *div(minimum(probability_E_excitation.Argument), 10):5:maximum(probability_E_excitation.Argument))
    Process_plot(plot_P_E, "probability_E_excitation", fissionant_nucleus_identifier)
end
if neutron_spectrum == "YES"
    plot_neutron_spectrum = Plot_data(n_E.E, Ratio_to_Maxwellian, "", :red)
    Modify_plot(
        plot_neutron_spectrum, "E [MeV]", "Nuetron spectrum 1/MeV", (first(n_E.E), last(n_E.E)),
        :identity, (minimum(n_E.Value), maximum(n_E.Value)), :log10, "Neutron spectrum, linear scale"
    )
end