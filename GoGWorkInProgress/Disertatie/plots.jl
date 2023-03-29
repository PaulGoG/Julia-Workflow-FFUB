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
function Plot_log10_yticks(plt::Plots.Plot)
    ticks = [10.0^i for i in -10:10]
    yticks!(
        plt,
        ticks
    )
end
function Plot_P_q_xticks(plt::Plots.Plot, maxval)
    xticks!(
        plt,
        0.0:0.25:maxval
    )
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
    Modify_plot(
        plot_y_A_log, "A", "Yield %", 
        (minimum(y_A.Argument), maximum(y_A.Argument)), :identity, 
        (minimum(y_A.Value), maximum(y_A.Value)*2.5), :log10, "Y(A) in log10 scale"
    )
    xticks!(plot_y_A_log, 10 *div(minimum(y_A.Argument), 10):5:maximum(y_A.Argument))
    Plot_log10_yticks(plot_y_A_log)
    Process_plot(plot_y_A_log, "y_A_log", fissionant_nucleus_identifier)

    plot_y_A_lin = Scatter_data(y_A.Argument, y_A.Value, "", :red, 5, :circle)
    Plot_data(plot_y_A_lin, y_A.Argument, y_A.Value, "", :red)
    Modify_plot(plot_y_A_lin)
    Modify_plot(
        plot_y_A_lin, "A", "Yield %", 
        (minimum(y_A.Argument), maximum(y_A.Argument)), :identity, 
        (0.0, maximum(y_A.Value)*1.1), :identity, "Y(A) in linear scale"
    )
    avg_A_L = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .<= A_H_min])
    avg_A_H = Average_yield_argument(y_A, y_A.Argument[y_A.Argument .>= A_H_min])
    Plot_textbox(plot_y_A_lin, A_H_min, maximum(y_A.Value)*1.05 - 0.5, latexstring("\$\\mathrm{<A>_L}\$ = $(round(avg_A_L[1], digits = 1))"))
    Plot_textbox(plot_y_A_lin, A_H_min, maximum(y_A.Value)*1.05, latexstring("\$\\mathrm{<A>_H}\$ = $(round(avg_A_H[1], digits = 1))"))
    xticks!(plot_y_A_lin, 10 *div(minimum(y_A.Argument), 10):5:maximum(y_A.Argument))
    Process_plot(plot_y_A_lin, "y_A_lin", fissionant_nucleus_identifier)

    plot_y_Z = Scatter_data(y_Z.Argument, y_Z.Value, "", :red, 5, :circle)
    Plot_data(plot_y_Z, y_Z.Argument, y_Z.Value, "", :red)
    Modify_plot(plot_y_Z)
    Modify_plot(
        plot_y_Z, "Z", "Yield %", 
        (minimum(y_Z.Argument), maximum(y_Z.Argument)), :identity, 
        (0.0, maximum(y_Z.Value)*1.1), :identity, ""
    )
    δₑₒ = (sum(y_Z.Value[iseven.(y_Z.Argument)]) - sum(y_Z.Value[isodd.(y_Z.Argument)]))/sum(y_Z.Value)
    Plot_textbox(plot_y_Z, Z₀/2, maximum(y_Z.Value)*1.05, "δₑₒ = $(round(δₑₒ *100, digits = 2)) %")
    xticks!(plot_y_Z, 10 *div(minimum(y_Z.Argument), 10):2:maximum(y_Z.Argument))
    Process_plot(plot_y_Z, "y_Z", fissionant_nucleus_identifier)

    plot_y_N = Scatter_data(y_N.Argument, y_N.Value, "", :red, 4, :circle)
    Modify_plot(plot_y_N)
    Modify_plot(
        plot_y_N, "N", "Yield %", 
        (minimum(y_N.Argument), maximum(y_N.Argument)), :identity, 
        (0.0, maximum(y_N.Value)*1.1), :identity, ""
    )
    xticks!(plot_y_N, 10 *div(minimum(y_N.Argument), 10):2:maximum(y_N.Argument))
    Process_plot(plot_y_N, "y_N", fissionant_nucleus_identifier)

    avg_TKE = Average_yield_argument(y_TKE, y_TKE.Argument)
    plot_y_TKE = Plot_data(y_TKE.Argument, y_TKE.Value, y_TKE.σ, "<TKE> = $(round(avg_TKE[1], digits=2)) MeV", :red)
    Modify_plot(plot_y_TKE)
    Modify_plot(
        plot_y_TKE, "TKE", "Yield %", 
        (minimum(y_TKE.Argument), maximum(y_TKE.Argument)), :identity, 
        (0.0, maximum(y_TKE.Value)*1.1), :identity, ""
    )
    xticks!(plot_y_TKE, 10 *div(minimum(y_TKE.Argument), 10):5:maximum(y_TKE.Argument))
    Process_plot(plot_y_TKE, "y_TKE", fissionant_nucleus_identifier)

    plot_tke_AH = Scatter_data(tke_AH.Argument, tke_AH.Value, "", :red, 5, :circle)
    Plot_data(plot_tke_AH, tke_AH.Argument, tke_AH.Value, "", :red)
    Modify_plot(plot_tke_AH)
    Modify_plot(
        plot_tke_AH, L"\mathrm{A_H}", "TKE [MeV]", 
        (minimum(tke_AH.Argument), maximum(tke_AH.Argument)), :identity, 
        (minimum(tke_AH.Value)*0.9, maximum(tke_AH.Value)*1.1), :identity, ""
    )
    xticks!(plot_tke_AH, 10 *div(minimum(tke_AH.Argument), 10):5:maximum(tke_AH.Argument))
    Process_plot(plot_tke_AH, "tke_AH", fissionant_nucleus_identifier)

    plot_ke_A = Scatter_data(ke_A.Argument, ke_A.Value, "", :red, 5, :circle)
    Plot_data(plot_ke_A, ke_A.Argument, ke_A.Value, "", :red)
    Modify_plot(plot_ke_A)
    Modify_plot(
        plot_ke_A, "A", "KE [MeV]", 
        (minimum(ke_A.Argument), maximum(ke_A.Argument)), :identity, 
        (minimum(ke_A.Value)*0.9, maximum(ke_A.Value)*1.1), :identity, ""
    )
    xticks!(plot_ke_A, 10 *div(minimum(ke_A.Argument), 10):5:maximum(ke_A.Argument))
    Process_plot(plot_ke_A, "ke_A", fissionant_nucleus_identifier)
end
if secondary_output_ν == "YES"
    #=
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
    =#
    gr(size = plots_resolution, dpi=300)

    avg_ν_L = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_L_range)
    avg_ν_H = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_H_range)
    plot_ν_A = Scatter_data(ν_A.Argument, ν_A.Value, "", :red, 5, :circle)
    Plot_data(plot_ν_A, ν_A.Argument, ν_A.Value, "", :red)
    Modify_plot(plot_ν_A)
    Modify_plot(
        plot_ν_A, "A", L"\nu", 
        (minimum(ν_A.Argument), maximum(ν_A.Argument)), :identity, 
        (0.0, maximum(ν_A.Value)*1.1), :identity, "Neutron multiplicity for all fragments"
    )
    Plot_textbox(plot_ν_A, ν_A.Argument[15], maximum(ν_A.Value)*1.05, latexstring("\$\\mathrm{<\\nu>_L}\$ = $(round(avg_ν_L, digits = 2))"))
    Plot_textbox(plot_ν_A, ν_A.Argument[15], maximum(ν_A.Value)*1.05 - 0.25, latexstring("\$\\mathrm{<\\nu>_H}\$ = $(round(avg_ν_H, digits = 2))"))
    xticks!(plot_ν_A, 10 *div(minimum(ν_A.Argument), 10):5:maximum(ν_A.Argument))
    yticks!(plot_ν_A, 0.0:0.5:maximum(ν_A.Value)*1.05)
    Process_plot(plot_ν_A, "nu_A", fissionant_nucleus_identifier)

    avg_ν_Pair = Average_value(ν_AH_Pair, y_A, A_H_range)
    plot_ν_Pair_AH = Scatter_data(ν_AH_Pair.Argument, ν_AH_Pair.Value, "", :red, 5, :circle)
    Plot_data(plot_ν_Pair_AH, ν_AH_Pair.Argument, ν_AH_Pair.Value, latexstring("\$\\mathrm{<\\nu_{pair}>}\$ = $(round(avg_ν_Pair, digits = 2))"), :red)
    Modify_plot(plot_ν_Pair_AH)
    Modify_plot(
        plot_ν_Pair_AH, L"\mathrm{A_H}", latexstring("\$\\mathrm{\\nu_{pair}}\$"), 
        (minimum(ν_AH_Pair.Argument), maximum(ν_AH_Pair.Argument)), :identity, 
        (minimum(ν_AH_Pair.Value)*0.8, maximum(ν_AH_Pair.Value)*1.2), :identity, "Neutron multiplicity for all fragment pairs"
    )
    xticks!(plot_ν_Pair_AH, 10 *div(minimum(ν_AH_Pair.Argument), 10):5:maximum(ν_AH_Pair.Argument))
    Process_plot(plot_ν_Pair_AH, "nuPair_AH", fissionant_nucleus_identifier)

    Ratio_ν = ν_A.Value[ν_A.Argument .>= A_H_min] ./ ν_AH_Pair.Value
    plot_Ratio_ν_AH = Scatter_data(ν_AH_Pair.Argument, Ratio_ν, "", :red, 5, :circle)
    Plot_data(plot_Ratio_ν_AH, ν_AH_Pair.Argument, Ratio_ν, "", :red)
    Modify_plot(plot_Ratio_ν_AH)
    Modify_plot(
        plot_Ratio_ν_AH, L"\mathrm{A_H}", L"\mathrm{\nu_H/\nu_{pair}}", 
        (minimum(ν_AH_Pair.Argument), maximum(ν_AH_Pair.Argument)), :identity, 
        (0.0, 1.0), :identity, "Neutron multiplicity ratio over HF domain"
    )
    hline!(plot_Ratio_ν_AH, [0.5], linestyle = :dashdot, color = :black, label = "")
    xticks!(plot_Ratio_ν_AH, 10 *div(minimum(ν_AH_Pair.Argument), 10):5:maximum(ν_AH_Pair.Argument))
    Process_plot(plot_Ratio_ν_AH, "Ratio_nuH_nuPair", fissionant_nucleus_identifier)

    avg_ν = Average_value(ν_A_Z_TKE, y_A_Z_TKE, A_range)
    plot_P_ν = Scatter_data(probability_ν.Argument, probability_ν.Value, "", :red, 5, :circle)
    Plot_data(plot_P_ν, probability_ν.Argument, probability_ν.Value, latexstring("\$<\\nu>\$ = $(round(avg_ν, digits=2))"), :red)
    Modify_plot(plot_P_ν)
    Modify_plot(
        plot_P_ν, L"\nu", "Probability %", 
        (minimum(probability_ν.Argument), maximum(probability_ν.Argument)), :identity, 
        (minimum(probability_ν.Value)*0.9, maximum(probability_ν.Value)*1.1), :identity, latexstring("Probability of occurance \$\\mathrm{P(\\nu)}\$")
    )
    Process_plot(plot_P_ν, "P_nu", fissionant_nucleus_identifier)

    if secondary_output_Ap == "YES"
        #=
        plotlyjs(size = (16, 9) .* 90, dpi=300)
        plot_heatmap_y_Ap_Z = Plot_heatmap(
            DataFrame(x = y_Ap_Z.A, y = y_Ap_Z.Z, z = y_Ap_Z.Value),
            10, 10, 2, 2,
            "Y(Aₚ,Z)", :turbo
        )
        Modify_plot(plot_heatmap_y_Ap_Z, "Aₚ", "Z", "")
        Modify_plot(plot_heatmap_y_Ap_Z)
        display(plot_heatmap_y_Ap_Z)
        =#
        gr(size = plots_resolution, dpi=300)

        avg_Ap_L = Average_yield_argument(y_Ap, y_Ap.Argument[y_Ap.Argument .< Ap_H_min])
        avg_Ap_H = Average_yield_argument(y_Ap, y_Ap.Argument[y_Ap.Argument .>= Ap_H_min])
        plot_y_Ap_lin = Scatter_data(y_Ap.Argument, y_Ap.Value, "", :red, 4, :circle)
        Plot_data(plot_y_Ap_lin, y_Ap.Argument, y_Ap.Value, "", :red)
        Modify_plot(plot_y_Ap_lin)
        Modify_plot(
            plot_y_Ap_lin, L"\mathrm{A_p}", "Yield %", 
            (minimum(y_Ap.Argument), maximum(y_Ap.Argument)), :identity, 
            (0.0, maximum(y_Ap.Value)*1.1), :identity, latexstring("\$\\mathrm{Y(A_p)}\$ in linear scale")
        )
        Plot_textbox(plot_y_Ap_lin, A_H_min, maximum(y_Ap.Value)*1.05 - 0.5, latexstring("\$\\mathrm{<A_p>_L}\$ = $(round(avg_Ap_L[1], digits = 1))"))
        Plot_textbox(plot_y_Ap_lin, A_H_min, maximum(y_Ap.Value)*1.05, latexstring("\$\\mathrm{<A_p>_H}\$ = $(round(avg_Ap_H[1], digits = 1))"))
        xticks!(plot_y_Ap_lin, 10 *div(minimum(y_Ap.Argument), 10):5:maximum(y_Ap.Argument))
        Process_plot(plot_y_Ap_lin, "y_Ap_lin", fissionant_nucleus_identifier)

        plot_y_Ap_comp_lin = Scatter_data(y_Ap.Argument, y_Ap.Value, L"\mathrm{Y(A_p)}", :red, 4, :circle)
        Plot_data(plot_y_Ap_comp_lin, y_Ap.Argument, y_Ap.Value, "", :red)
        Scatter_data(plot_y_Ap_comp_lin, y_A.Argument, y_A.Value, L"\mathrm{Y(A)}", :blue, 4, :xcross)
        Plot_data(plot_y_Ap_comp_lin, y_A.Argument, y_A.Value, "", :blue)
        Modify_plot(plot_y_Ap_comp_lin)
        Modify_plot(
            plot_y_Ap_comp_lin, L"\mathrm{A_p}", "Yield %", 
            (minimum(y_Ap.Argument), maximum(y_Ap.Argument)), :identity, 
            (0.0, maximum(y_Ap.Value)*1.1), :identity, latexstring("\$\\mathrm{Y(A)\\: and\\: Y(A_p)}\$ in linear scale")
        )
        xticks!(plot_y_Ap_lin, 10 *div(minimum(y_Ap.Argument), 10):5:maximum(y_Ap.Argument))
        Process_plot(plot_y_Ap_comp_lin, "y_Ap_comparisson_lin", fissionant_nucleus_identifier)
        
        plot_y_Ap_log = Scatter_data(y_Ap.Argument, y_Ap.Value, "", :red, 4, :circle)
        Plot_data(plot_y_Ap_log, y_Ap.Argument, y_Ap.Value, "", :red)
        Modify_plot(plot_y_Ap_log)
        Modify_plot(
            plot_y_Ap_log, L"\mathrm{A_p}", "Yield %", 
            (minimum(y_Ap.Argument), maximum(y_Ap.Argument)), :identity, 
            (minimum(y_Ap.Value), maximum(y_Ap.Value)*2), :log10, latexstring("\$\\mathrm{Y(A_p)}\$ in log10 scale")
        )
        xticks!(plot_y_Ap_log, 10 *div(minimum(y_Ap.Argument), 10):5:maximum(y_Ap.Argument))
        Plot_log10_yticks(plot_y_Ap_log)
        Process_plot(plot_y_Ap_log, "y_Ap_log", fissionant_nucleus_identifier)

        plot_y_Ap_comp_log = Scatter_data(y_Ap.Argument, y_Ap.Value, L"\mathrm{Y(A_p)}", :red, 4, :circle)
        Plot_data(plot_y_Ap_comp_log, y_Ap.Argument, y_Ap.Value, "", :red)
        Scatter_data(plot_y_Ap_comp_log, y_A.Argument, y_A.Value, L"\mathrm{Y(A)}", :blue, 4, :xcross)
        Plot_data(plot_y_Ap_comp_log, y_A.Argument, y_A.Value, "", :blue)
        Modify_plot(plot_y_Ap_comp_log)
        Modify_plot(
            plot_y_Ap_comp_log, L"\mathrm{A_p}", "Yield %", 
            (minimum(y_Ap.Argument), maximum(y_Ap.Argument)), :identity, 
            (minimum(y_Ap.Value), maximum(y_Ap.Value)*2), :log10, latexstring("\$\\mathrm{Y(A)\\: and\\: Y(A_p)}\$ in log10 scale")
        )
        xticks!(plot_y_Ap_comp_log, 10 *div(minimum(y_Ap.Argument), 10):5:maximum(y_Ap.Argument))
        Plot_log10_yticks(plot_y_Ap_comp_log)
        Plot_legend_attributes(plot_y_Ap_comp_log, :bottomright)
        Process_plot(plot_y_Ap_comp_log, "y_Ap_comparisson_log", fissionant_nucleus_identifier)

        plot_y_N_comp = Scatter_data(y_N.Argument, y_N.Value, L"\mathrm{Y(N)}", :blue, 4, :xcross)
        Scatter_data(plot_y_N_comp, y_Np.Argument, y_Np.Value, L"\mathrm{Y(N_p)}", :red, 4, :circle)
        Modify_plot(plot_y_N_comp)
        xticks!(plot_y_N_comp, 10 *div(minimum(y_Np.Argument), 10):2:maximum(y_Np.Argument))
        Modify_plot(
            plot_y_N_comp, "N", "Yield %", 
            (minimum(y_Np.Argument), maximum(y_Np.Argument)), :identity, 
            (0.0, maximum(y_Np.Value)*1.1), :identity, ""
        )
        Process_plot(plot_y_N_comp, "y_Np_comparisson", fissionant_nucleus_identifier)
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
    plot_P_T = Bar_data(probability_T.Argument, probability_T.Value, latexstring("\$\\mathrm{<T>}\$ = $(round(avg_T, digits = 2)) MeV"), :red, ΔT)
    Modify_plot(plot_P_T)
    Modify_plot(
        plot_P_T, "T [MeV]", "Probability %", 
        (0.0, maximum(probability_T.Argument)), :identity, 
        (0.0, maximum(probability_T.Value)*1.1), :identity, latexstring("Probability of occurance \$\\mathrm{P(T)}\$ for all fragments")
    )
    Plot_P_q_xticks(plot_P_T, maximum(probability_T.Argument))
    Process_plot(plot_P_T, "P_T", fissionant_nucleus_identifier)

    avg_T_L = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_L_range)
    plot_P_T_L = Bar_data(probability_T_L.Argument, probability_T_L.Value, latexstring("\$\\mathrm{<T>_L}\$ = $(round(avg_T_L, digits = 2)) MeV"), :red, ΔT)
    Modify_plot(plot_P_T_L)
    Modify_plot(
        plot_P_T_L, "T [MeV]", "Probability %", 
        (0.0, maximum(probability_T_L.Argument)), :identity, 
        (0.0, maximum(probability_T_L.Value)*1.1), :identity, latexstring("Probability of occurance \$\\mathrm{P(T_L)}\$ for light fragments")
    )
    Plot_P_q_xticks(plot_P_T_L, maximum(probability_T_L.Argument))
    Process_plot(plot_P_T_L, "P_T_LF", fissionant_nucleus_identifier)

    avg_T_H = Average_value(T_A_Z_TKE, y_A_Z_TKE, A_H_range)
    plot_P_T_H = Bar_data(probability_T_H.Argument, probability_T_H.Value, latexstring("\$\\mathrm{<T>_H}\$ = $(round(avg_T_H, digits = 2)) MeV"), :red, ΔT)
    Modify_plot(plot_P_T_H)
    Modify_plot(
        plot_P_T_H, "T [MeV]", "Probability %", 
        (0.0, maximum(probability_T_H.Argument)), :identity, 
        (0.0, maximum(probability_T_H.Value)*1.1), :identity, latexstring("Probability of occurance \$\\mathrm{P(T_H)}\$ for heavy fragments")
    )
    Plot_P_q_xticks(plot_P_T_H, maximum(probability_T_H.Argument))
    Process_plot(plot_P_T_H, "P_T_HF", fissionant_nucleus_identifier)
end
if secondary_output_avg_ε == "YES"
    gr(size = plots_resolution, dpi=300)

    avg_ε = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_range)
    plot_avgE_SCM = Bar_data(probability_avg_ε.Argument, probability_avg_ε.Value, latexstring("\$\\overline{<\\varepsilon>}\$ = $(round(avg_ε, digits = 2)) MeV"), :red, Δavg_ε)
    Modify_plot(plot_avgE_SCM)
    Modify_plot(
        plot_avgE_SCM, L"\mathrm{<\varepsilon>\: [MeV]}", "Probability %", 
        (0.0, maximum(probability_avg_ε.Argument)), :identity, 
        (0.0, maximum(probability_avg_ε.Value)*1.1), :identity, latexstring("Probability of occurance \$\\mathrm{P(<\\varepsilon>)}\$ for all fragments")
    )
    Plot_P_q_xticks(plot_avgE_SCM, maximum(probability_avg_ε.Argument))
    Process_plot(plot_avgE_SCM, "P_avgE_SCM", fissionant_nucleus_identifier)

    avg_ε_L = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_L_range)
    plot_avgE_L_SCM = Bar_data(probability_avg_ε_L.Argument, probability_avg_ε_L.Value, latexstring("\$\\overline{<\\varepsilon>_L}\$ = $(round(avg_ε_L, digits = 2)) MeV"), :red, Δavg_ε)
    Modify_plot(plot_avgE_L_SCM)
    Modify_plot(
        plot_avgE_L_SCM, L"\mathrm{<\varepsilon>\:[MeV]}", "Probability %", 
        (0.0, maximum(probability_avg_ε_L.Argument)), :identity, 
        (0.0, maximum(probability_avg_ε_L.Value)*1.1), :identity, latexstring("Probability of occurance \$\\mathrm{P(<\\varepsilon>_L)}\$ for light fragments")
    )
    Plot_P_q_xticks(plot_avgE_L_SCM, maximum(probability_avg_ε_L.Argument))
    Process_plot(plot_avgE_L_SCM, "P_avgE_LF_SCM", fissionant_nucleus_identifier)

    avg_ε_H = Average_value(avg_ε_A_Z_TKE, y_A_Z_TKE, A_H_range)
    plot_avgE_H_SCM = Bar_data(probability_avg_ε_H.Argument, probability_avg_ε_H.Value, latexstring("\$\\overline{<\\varepsilon>_H}\$ = $(round(avg_ε_H, digits = 2)) MeV"), :red, Δavg_ε)
    Modify_plot(plot_avgE_H_SCM)
    Modify_plot(
        plot_avgE_H_SCM, L"\mathrm{<\varepsilon>\: [MeV]}", "Probability %", 
        (0.0, maximum(probability_avg_ε_H.Argument)), :identity, 
        (0.0, maximum(probability_avg_ε_H.Value)*1.1), :identity, latexstring("Probability of occurance \$\\mathrm{P(<\\varepsilon>_H)}\$ for heavy fragments")
    )
    Plot_P_q_xticks(plot_avgE_H_SCM, maximum(probability_avg_ε_H.Argument))
    Process_plot(plot_avgE_H_SCM, "P_avgE_HF_SCM", fissionant_nucleus_identifier)
end
if secondary_output_TXE_Q == "YES"
    gr(size = plots_resolution, dpi=300)

    avg_Q = Average_value(Q_AH, y_A, A_H_range)
    plot_Q_AH = Scatter_data(Q_AH.Argument, Q_AH.Value, "", :red, 5, :circle)
    Plot_data(plot_Q_AH, Q_AH.Argument, Q_AH.Value, "<Q> = $(round(avg_Q, digits = 2)) MeV", :red)
    Modify_plot(plot_Q_AH)
    Modify_plot(
        plot_Q_AH, L"\mathrm{A_H}", "Q [MeV]", 
        (minimum(Q_AH.Argument), maximum(Q_AH.Argument)), :identity, 
        (minimum(Q_AH.Value)*0.9, maximum(Q_AH.Value)*1.1), :identity, ""
    )
    xticks!(plot_Q_AH, 10 *div(minimum(Q_AH.Argument), 10):5:maximum(Q_AH.Argument))
    Process_plot(plot_Q_AH, "Q_AH", fissionant_nucleus_identifier)

    avg_TXE = Average_value(txe_AH, y_A, A_H_range)
    plot_TXE_AH = Scatter_data(txe_AH.Argument, txe_AH.Value, "", :red, 5, :circle)
    Plot_data(plot_TXE_AH, txe_AH.Argument, txe_AH.Value, "<TXE> = $(round(avg_TXE, digits = 2)) MeV", :red)
    Modify_plot(plot_TXE_AH)
    Modify_plot(
        plot_TXE_AH, L"\mathrm{A_H}", "TXE [MeV]", 
        (minimum(txe_AH.Argument), maximum(txe_AH.Argument)), :identity, 
        (minimum(txe_AH.Value)*0.9, maximum(txe_AH.Value)*1.1), :identity, ""
    )
    xticks!(plot_TXE_AH, 10 *div(minimum(txe_AH.Argument), 10):5:maximum(txe_AH.Argument))
    Process_plot(plot_TXE_AH, "TXE_AH", fissionant_nucleus_identifier)
end
if secondary_output_E_excitation == "YES"
    gr(size = plots_resolution, dpi=300)

    avg_E_exi = Average_value(E_excitation_A, y_A, A_range)
    plot_P_E = Bar_data(probability_E_excitation.Argument, probability_E_excitation.Value, "<E*> = $(round(avg_E_exi, digits = 2)) MeV", :red, 1.0)
    Modify_plot(plot_P_E)
    Modify_plot(
        plot_P_E, "E* [MeV]", "Probability %", 
        (0.0, maximum(probability_E_excitation.Argument)/2), :identity, 
        (0.0, maximum(probability_E_excitation.Value)*1.1), :identity, "Probability of occurance P(E*) for initial fragments"
    )
    xticks!(plot_P_E, 0.0:2.5:maximum(probability_E_excitation.Argument))
    Process_plot(plot_P_E, "P_E_exi", fissionant_nucleus_identifier)
end
if neutron_spectrum == "YES"
    plot_nSpectrum = Plot_data(n_E.E, n_E.Value, "", :red)
    Modify_plot(
        plot_nSpectrum, "E [MeV]", "Neutron spectrum 1/MeV", (0.0, last(n_E.E)),
        :identity, (0.0, maximum(n_E.Value)*1.1), :identity, ""
    )
    Modify_plot(plot_nSpectrum)
    xticks!(plot_nSpectrum, 0.0:2.5:last(n_E.E))
    Process_plot(plot_nSpectrum, "nSpectrum_linscale", fissionant_nucleus_identifier)
    Modify_plot(
        plot_nSpectrum, "E [MeV]", "Neutron spectrum 1/MeV", (0.0, last(n_E.E)),
        :identity, (minimum(n_E.Value) *5e-1, maximum(n_E.Value) *10), :log10, ""
    )
    Plot_log10_yticks(plot_nSpectrum)
    Process_plot(plot_nSpectrum, "nSpectrum_logscale", fissionant_nucleus_identifier)

    plot_nSpectrum_ratio_Maxwellian = Plot_data(n_E.E, Ratio_to_Maxwellian, "", :red)
    Modify_plot(
        plot_nSpectrum_ratio_Maxwellian, "E [MeV]", "Neutron spectrum as ratio to Maxwellian", (first(n_E.E), 10.0),
        :log10, (0.5, 1.5), :identity, ""
    )
    hline!(plot_nSpectrum_ratio_Maxwellian, [1.0], linestyle = :dashdot, color = :black, label = latexstring("\$\\mathrm{T_M}\$ = $(round(T_M_eq, digits = 2)) MeV"))
    Modify_plot(plot_nSpectrum_ratio_Maxwellian)
    xticks!(plot_nSpectrum_ratio_Maxwellian, [10.0^i for i in -10:10])
    Process_plot(plot_nSpectrum_ratio_Maxwellian, "nSpectrum_ratio_Maxwellian", fissionant_nucleus_identifier)
end