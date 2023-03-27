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

    plot_y_A_log = Scatter_data(y_A.Argument, y_A.Value, "", :red, 5, :circle)
    Plot_data(plot_y_A_log, y_A.Argument, y_A.Value, "", :red)
    Modify_plot(plot_y_A_log)
    xticks!(plot_y_A_log, 10 *div(minimum(y_A.Argument), 10):5:maximum(y_A.Argument))
    Modify_plot(
        plot_y_A_log, "A", "Yield %", 
        (minimum(y_A.Argument), maximum(y_A.Argument)), :identity, 
        (minimum(y_A.Value)*5e-1, maximum(y_A.Value)*2.5), :log10, ""
    )
    Process_plot(plot_y_A_log, "y_A_log", fissionant_nucleus_identifier)

    plot_y_A_lin = Scatter_data(y_A.Argument, y_A.Value, "", :red, 5, :circle)
    Plot_data(plot_y_A_lin, y_A.Argument, y_A.Value, "", :red)
    Modify_plot(plot_y_A_lin)
    xticks!(plot_y_A_lin, 10 *div(minimum(y_A.Argument), 10):5:maximum(y_A.Argument))
    Modify_plot(
        plot_y_A_lin, "A", "Yield %", 
        (minimum(y_A.Argument), maximum(y_A.Argument)), :identity, 
        (minimum(y_A.Value)*0.9, maximum(y_A.Value)*1.1), :identity, ""
    )
    avg_A_L = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .<= A_H_min])
    avg_A_H = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .>= A_H_min])
    Plot_textbox(plot_y_A_lin, A_H_min, maximum(y_A.Value)*1.05, latexstring("\$\\mathrm{<A_L>}\$ = $(round(avg_A_L[1], digits = 1))"))
    Plot_textbox(plot_y_A_lin, A_H_min, maximum(y_A.Value)*1.05 - 0.5, latexstring("\$\\mathrm{<A_H>}\$ = $(round(avg_A_H[1], digits = 1))"))
    Process_plot(plot_y_A_lin, "y_A_lin", fissionant_nucleus_identifier)

    plot_y_Z = Scatter_data(y_Z.Argument, y_Z.Value, "", :red, 5, :circle)
    Plot_data(plot_y_Z, y_Z.Argument, y_Z.Value, "", :red)
    Modify_plot(plot_y_Z)
    xticks!(plot_y_Z, 10 *div(minimum(y_Z.Argument), 10):2:maximum(y_Z.Argument))
    Modify_plot(
        plot_y_Z, "Z", "Yield %", 
        (minimum(y_Z.Argument), maximum(y_Z.Argument)), :identity, 
        (minimum(y_Z.Value)*0.9, maximum(y_Z.Value)*1.1), :identity, ""
    )
    Process_plot(plot_y_Z, "y_Z", fissionant_nucleus_identifier)

    plot_y_N = Scatter_data(y_N.Argument, y_N.Value, "", :red, 5, :circle)
    Plot_data(plot_y_N, y_N.Argument, y_N.Value, "", :red)
    Modify_plot(plot_y_N)
    xticks!(plot_y_N, 10 *div(minimum(y_N.Argument), 10):2:maximum(y_N.Argument))
    Modify_plot(
        plot_y_N, "N", "Yield %", 
        (minimum(y_N.Argument), maximum(y_N.Argument)), :identity, 
        (minimum(y_N.Value)*0.9, maximum(y_N.Value)*1.1), :identity, ""
    )
    Process_plot(plot_y_N, "y_N", fissionant_nucleus_identifier)

    plot_y_TKE = Scatter_data(y_TKE.Argument, y_TKE.Value, "", :red, 5, :circle)
    Plot_data(plot_y_TKE, y_TKE.Argument, y_TKE.Value, "", :red)
    Modify_plot(plot_y_TKE)
    xticks!(plot_y_TKE, 10 *div(minimum(y_TKE.Argument), 10):5:maximum(y_TKE.Argument))
    Modify_plot(
        plot_y_TKE, "TKE", "Yield %", 
        (minimum(y_TKE.Argument), maximum(y_TKE.Argument)), :identity, 
        (minimum(y_TKE.Value)*0.9, maximum(y_TKE.Value)*1.1), :identity, ""
    )
    Process_plot(plot_y_TKE, "y_TKE", fissionant_nucleus_identifier)

    plot_tke_AH = Scatter_data(tke_AH.Argument, tke_AH.Value, "", :red, 5, :circle)
    Plot_data(plot_tke_AH, tke_AH.Argument, tke_AH.Value, "", :red)
    Modify_plot(plot_tke_AH)
    xticks!(plot_tke_AH, 10 *div(minimum(tke_AH.Argument), 10):5:maximum(tke_AH.Argument))
    Modify_plot(
        plot_tke_AH, L"\mathrm{A_H}", "TKE [MeV]", 
        (minimum(tke_AH.Argument), maximum(tke_AH.Argument)), :identity, 
        (minimum(tke_AH.Value)*0.9, maximum(tke_AH.Value)*1.1), :identity, ""
    )
    Process_plot(plot_tke_AH, "tke_AH", fissionant_nucleus_identifier)

    plot_ke_A = Scatter_data(ke_A.Argument, ke_A.Value, "", :red, 5, :circle)
    Plot_data(plot_ke_A, ke_A.Argument, ke_A.Value, "", :red)
    Modify_plot(plot_ke_A)
    xticks!(plot_ke_A, 10 *div(minimum(ke_A.Argument), 10):5:maximum(ke_A.Argument))
    Modify_plot(
        plot_ke_A, "A", "KE [MeV]", 
        (minimum(ke_A.Argument), maximum(ke_A.Argument)), :identity, 
        (minimum(ke_A.Value)*0.9, maximum(ke_A.Value)*1.1), :identity, ""
    )
    Process_plot(plot_ke_A, "ke_A", fissionant_nucleus_identifier)
end
if secondary_output_ν == "YES"
    plotlyjs(size = (16, 9) .* 90, dpi=600)

    plot_surface_ν_A_TKE = Plot_surface(
        DataFrame(x = ν_A_TKE.A, y = ν_A_TKE.TKE, z = ν_A_TKE.Value),
        10, 10, 20, Int(round(10/TKE_step)), (120, 0), 
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
        (minimum(ν_A.Value)*0.9, maximum(ν_A.Value)*1.1), :identity, ""
    )
    Plot_textbox(plot_ν_A, ν_A.Argument[15], maximum(ν_A.Value)*1.05, latexstring("\$\\mathrm{<\\nu_L>}\$ = $(round(avg_ν_L, digits = 2))"))
    Plot_textbox(plot_ν_A, ν_A.Argument[15], maximum(ν_A.Value)*1.05 - 0.25, latexstring("\$\\mathrm{<\\nu_H>}\$ = $(round(avg_ν_H, digits = 2))"))
    xticks!(plot_ν_A, 10 *div(minimum(ν_A.Argument), 10):5:maximum(ν_A.Argument))
    Process_plot(plot_ν_A, "nu_A", fissionant_nucleus_identifier)

    avg_ν_Pair = Average_value(ν_AH_Pair, y_A, A_H_range)
    plot_ν_Pair_AH = Scatter_data(ν_AH_Pair.Argument, ν_AH_Pair.Value, "", :red, 5, :circle)
    Plot_data(plot_ν_Pair_AH, ν_AH_Pair.Argument, ν_AH_Pair.Value, "", :red)
    Modify_plot(plot_ν_Pair_AH)
    Modify_plot(
        plot_ν_Pair_AH, L"\mathrm{A_H}", latexstring("\$\\mathrm{\\nu_{pair}}\$"), 
        (minimum(ν_AH_Pair.Argument), maximum(ν_AH_Pair.Argument)), :identity, 
        (minimum(ν_AH_Pair.Value)*0.9, maximum(ν_AH_Pair.Value)*1.1), :identity, ""
    )
    Plot_textbox(plot_ν_Pair_AH, ν_AH_Pair.Argument[10], maximum(ν_AH_Pair.Value)*1.1 - 0.25, latexstring("\$\\mathrm{<\\nu_{pair}>}\$ = $(round(avg_ν_Pair, digits = 2))"))
    xticks!(plot_ν_Pair_AH, 10 *div(minimum(ν_AH_Pair.Argument), 10):5:maximum(ν_AH_Pair.Argument))
    Process_plot(plot_ν_Pair_AH, "nuPair_AH", fissionant_nucleus_identifier)

    Ratio_ν = ν_A.Value[ν_A.Argument .>= A_H_min] ./ ν_AH_Pair.Value
    plot_Ratio_ν_AH = Scatter_data(ν_AH_Pair.Argument, Ratio_ν, "", :red, 5, :circle)
    Plot_data(plot_Ratio_ν_AH, ν_AH_Pair.Argument, Ratio_ν, "", :red)
    Modify_plot(plot_Ratio_ν_AH)
    Modify_plot(
        plot_Ratio_ν_AH, L"\mathrm{A_H}", L"\nu_H/(\nu_L + \nu_H)", 
        (minimum(ν_AH_Pair.Argument), maximum(ν_AH_Pair.Argument)), :identity, 
        (0.0, 1.0), :identity, ""
    )
    hline!(plot_Ratio_ν_AH, [0.5], linestyle = :dashdot, color = :black, label = "")
    xticks!(plot_Ratio_ν_AH, 10 *div(minimum(ν_AH_Pair.Argument), 10):5:maximum(ν_AH_Pair.Argument))
    Process_plot(plot_Ratio_ν_AH, "Ratio_nuH_nuPair", fissionant_nucleus_identifier)

    plot_P_ν = Scatter_data(probability_ν.Argument, probability_ν.Value, "", :red, 5, :circle)
    Plot_data(plot_P_ν, probability_ν.Argument, probability_ν.Value, "", :red)
    Modify_plot(plot_P_ν)
    Modify_plot(
        plot_P_ν, "ν", "Probability %", 
        (minimum(probability_ν.Argument), maximum(probability_ν.Argument)), :identity, 
        (minimum(probability_ν.Value)*0.9, maximum(probability_ν.Value)*1.1), :identity, ""
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
            (minimum(y_Ap.Value)*0.9, maximum(y_Ap.Value)*1.1), :identity, ""
        )
        Plot_textbox(plot_y_Ap, A_H_min, maximum(y_Ap.Value)*1.05, latexstring("\$\\mathrm{<Ap_L>}\$ = $(round(avg_Ap_L[1], digits = 1))"))
        Plot_textbox(plot_y_Ap, A_H_min, maximum(y_Ap.Value)*1.05 - 0.5, latexstring("\$\\mathrm{<Ap_H>}\$ = $(round(avg_Ap_H[1], digits = 1))"))
        xticks!(plot_y_Ap, 10 *div(minimum(y_Ap.Argument), 10):5:maximum(y_Ap.Argument))
        Process_plot(plot_y_Ap, "y_Ap_lin", fissionant_nucleus_identifier)
        
        #MODIFY LIN VS LOG AND ADD COMPARISSONS W pre-neutron YIELD!
        Modify_plot(
            plot_y_Ap, "Aₚ", "Yield %", 
            (minimum(y_Ap.Argument), maximum(y_Ap.Argument)), :identity, 
            (minimum(y_Ap.Value)*5e-1, maximum(y_Ap.Value)*2), :log10, ""
        )
        Process_plot(plot_y_Ap, "y_Ap_log", fissionant_nucleus_identifier)
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

    avg_T = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_range)
    plot_P_T = Bar_data(probability_T.Argument, probability_T.Value, "", :red, ΔT)
    Modify_plot(plot_P_T)
    Modify_plot(
        plot_P_T, "T [MeV]", "Probability %", 
        (0.0, maximum(probability_T.Argument)), :identity, 
        (0.0, maximum(probability_T.Value)*1.1), :identity, ""
    )
    Plot_textbox(plot_P_T, 0.25, maximum(probability_T.Value)*1.07, "<T> = $(round(avg_T, digits = 2)) MeV")
    Process_plot(plot_P_T, "P_T", fissionant_nucleus_identifier)
end
if secondary_output_avg_ε == "YES"
    gr(size = plots_resolution, dpi=300)

    avg_ε = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_range)
    plot_avgE_SCM = Bar_data(probability_avg_ε.Argument, probability_avg_ε.Value, "", :red, Δavg_ε)
    Modify_plot(plot_avgE_SCM)
    Modify_plot(
        plot_avgE_SCM, "<ε> [MeV]", "Probability %", 
        (0.0, maximum(probability_avg_ε.Argument)), :identity, 
        (0.0, maximum(probability_avg_ε.Value)*1.1), :identity, ""
    )
    Plot_textbox(plot_avgE_SCM, 0.6, maximum(probability_avg_ε.Value)*1.075, "<ε>_avg = $(round(avg_ε, digits = 2)) MeV")
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
        (minimum(Q_AH.Value)*0.9, maximum(Q_AH.Value)*1.1), :identity, ""
    )
    Plot_textbox(plot_Q_AH, Q_AH.Argument[10], maximum(Q_AH.Value)*1.1 - 5, "<Q> = $(round(avg_Q, digits = 2)) MeV")
    xticks!(plot_Q_AH, 10 *div(minimum(Q_AH.Argument), 10):5:maximum(Q_AH.Argument))
    Process_plot(plot_Q_AH, "Q_AH", fissionant_nucleus_identifier)

    avg_TXE = Average_value(txe_AH, y_A, A_H_range)
    plot_TXE_AH = Scatter_data(txe_AH.Argument, txe_AH.Value, "", :red, 5, :circle)
    Plot_data(plot_TXE_AH, txe_AH.Argument, txe_AH.Value, "", :red)
    Modify_plot(plot_TXE_AH)
    Modify_plot(
        plot_TXE_AH, L"\mathrm{A_H}", "TXE [MeV]", 
        (minimum(txe_AH.Argument), maximum(txe_AH.Argument)), :identity, 
        (minimum(txe_AH.Value)*0.9, maximum(txe_AH.Value)*1.1), :identity, ""
    )
    Plot_textbox(plot_TXE_AH, txe_AH.Argument[10], maximum(txe_AH.Value)*1.075, "<TXE> = $(round(avg_TXE, digits = 2)) MeV")
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
        (0.0, maximum(probability_E_excitation.Argument)/2), :identity, 
        (0.0, maximum(probability_E_excitation.Value)*1.1), :identity, ""
    )
    Plot_textbox(plot_P_E, maximum(probability_E_excitation.Argument)/2 *0.85, maximum(probability_E_excitation.Value)*1.05, "<E*> = $(round(avg_E_exi, digits = 2)) MeV")
    xticks!(plot_P_E, 10 *div(minimum(probability_E_excitation.Argument), 10):5:maximum(probability_E_excitation.Argument))
    Process_plot(plot_P_E, "probability_E_excitation", fissionant_nucleus_identifier)
end
if neutron_spectrum == "YES"
    plot_nSpectrum = Plot_data(n_E.E, n_E.Value, "", :red)
    Modify_plot(
        plot_nSpectrum, "E [MeV]", "Neutron spectrum 1/MeV", (first(n_E.E), last(n_E.E)),
        :identity, (minimum(n_E.Value)*0.9, maximum(n_E.Value)*1.1), :identity, ""
    )
    Modify_plot(plot_nSpectrum)
    Process_plot(plot_nSpectrum, "nSpectrum_linscale", fissionant_nucleus_identifier)
    Modify_plot(
        plot_nSpectrum, "E [MeV]", "Neutron spectrum 1/MeV", (first(n_E.E), last(n_E.E)),
        :identity, (minimum(n_E.Value)*0.5, maximum(n_E.Value)*5), :log10, ""
    )
    Process_plot(plot_nSpectrum, "nSpectrum_logscale", fissionant_nucleus_identifier)

    plot_nSpectrum_ratio_Maxwellian = Plot_data(n_E.E, Ratio_to_Maxwellian, "", :red)
    Modify_plot(
        plot_nSpectrum_ratio_Maxwellian, "E [MeV]", "Neutron spectrum as ratio to Maxwellian", (first(n_E.E), 10.0),
        :log10, (0.5, 1.5), :identity, ""
    )
    hline!(plot_nSpectrum_ratio_Maxwellian, [1.0], linestyle = :dashdot, color = :black, label = "")
    Modify_plot(plot_nSpectrum_ratio_Maxwellian)
    Plot_textbox(plot_nSpectrum_ratio_Maxwellian, 1.0, 1.45, latexstring("\$\\mathrm{T_M}\$ = $(round(T_M_eq, digits = 2)) MeV"))
    Process_plot(plot_nSpectrum_ratio_Maxwellian, "nSpectrum_ratio_Maxwellian", fissionant_nucleus_identifier)
end