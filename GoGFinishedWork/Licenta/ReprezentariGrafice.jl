#= 
Sectiunea de cod unde se controleaza detaliile 
legate de generarea reprezentarilor grafice
=#
function Reprezinta_Suprafata(x, y, P, xinf, zheader, titlu)
    plt = surface(x, y, transpose(P),  
    xlabel = "x (km)", ylabel = "y (km)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    camera = (60,20), 
    c = :nipy_spectral, 
    yformatter = y->string(Int(y/1000)),
    xformatter = x->string(Int(x/1000)),
    zformatter = :scientific)
    display(plt)
    savefig("Reprezentari_Grafice\\SurfacePlot_$(titlu)_1.png")

    plt = surface(x, y, transpose(P),  
    xlabel = "x (km)", ylabel = "y (km)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    camera = (60,20), 
    c = :turbo, 
    yformatter = y->string(Int(y/1000)),
    xformatter = x->string(Int(x/1000)),
    zformatter = :scientific)
    display(plt)
    savefig("Reprezentari_Grafice\\SurfacePlot_$(titlu)_2.png")
end
function Reprezinta_Gradient(x, y, P, xinf, zheader, titlu)
    plt = heatmap(x, y, transpose(P),  
    xlabel = "x (km)", ylabel = "y (km)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    c = :turbo, 
    yformatter = y->string(Int(y/1000)),
    xformatter = x->string(Int(x/1000)),
    zformatter = :scientific,
    colorbar_title = zheader)
    display(plt)
    savefig("Reprezentari_Grafice\\HeatmapPlot_$(titlu)_1.png")

    plt = heatmap(x, y, transpose(P),  
    xlabel = "x (km)", ylabel = "y (km)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    c = :thermal, 
    yformatter = y->string(Int(y/1000)),
    xformatter = x->string(Int(x/1000)),
    zformatter = :scientific,
    colorbar_title = zheader)
    display(plt)
    savefig("Reprezentari_Grafice\\HeatmapPlot_$(titlu)_2.png")
end
function Reprezinta_Contur(x, y, P, xinf, zheader, titlu)
    plt = contourf(x, y, transpose(P),  
    xlabel = "x (km)", ylabel = "y (km)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    c = :turbo, 
    yformatter = y->string(Int(y/1000)),
    xformatter = x->string(Int(x/1000)),
    zformatter = :scientific,
    colorbar_title = zheader)
    display(plt)
    savefig("Reprezentari_Grafice\\ContourPlot_$(titlu)_1.png")

    plt = contourf(x, y, transpose(P),  
    xlabel = "x (km)", ylabel = "y (km)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    c = :thermal, 
    yformatter = y->string(Int(y/1000)),
    xformatter = x->string(Int(x/1000)),
    zformatter = :scientific,
    colorbar_title = zheader)
    display(plt)
    savefig("Reprezentari_Grafice\\ContourPlot_$(titlu)_2.png")
end
function Graficul_Resuspensiei(t_zile)
    y = [Coeficient_Resuspensie(i) for i in 1:t_zile]
    x = collect(1:1:t_zile)
    plt = plot(x, y,
    xlabel = "timp (zile)", 
    ylabel = "K (1/m)",
    color = "orange", 
    grid = false,
    plot_title = "K(t)",
    framestyle = :box,
    legend = false)
    display(plt)
    savefig("Reprezentari_Grafice\\CoeficientResuspensie.png")
end