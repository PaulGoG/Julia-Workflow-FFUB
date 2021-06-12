#= 
Sectiunea de cod unde se controleaza detaliile 
legate de generarea reprezentarilor grafice
=#
function Reprezinta_Suprafata(x, y, P, xinf, zheader, titlu)
    plt = surface(x, y, transpose(P),  
    xlabel = "x (m)", ylabel = "y (m)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    dpi = 1200, camera = (60,20), 
    c = :turbo)
    display(plt)
    savefig("Reprezentari_Grafice\\SurfacePlot_$(titlu).png")
end
function Reprezinta_Gradient(x, y, P, xinf, zheader, titlu)
    plt = heatmap(x, y, transpose(P),  
    xlabel = "x (m)", ylabel = "y (m)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    dpi = 1200)
    display(plt)
    savefig("Reprezentari_Grafice\\HeatmapPlot_$(titlu).png")
end
function Reprezinta_Contur(x, y, P, xinf, zheader, titlu)
    plt = contourf(x, y, transpose(P),  
    xlabel = "x (m)", ylabel = "y (m)", zlabel = zheader,
    xlim = (xinf, x[length(x)]),
    dpi = 1200)
    display(plt)
    savefig("Reprezentari_Grafice\\ContourPlot_$(titlu).png")
end