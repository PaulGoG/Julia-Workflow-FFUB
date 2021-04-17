# Plotarea unei matrici NxN
function Reprezinta_Suprafata(x, y, P)
    surface(x, y, transpose(P),  
    xlabel = "x", ylabel = "y",
    xlim = (0, x[length(x)]))
    savefig("Reprezentari_Grafice\\SurfacePlot.png")
end
function Reprezinta_Gradient(x, y, P)
    heatmap(x, y, transpose(P),  
    xlabel = "x", ylabel = "y")
    #savefig("Reprezentari_Grafice\\Denumire.png")
end
function Reprezinta_Contur(x, y, P)
    contourf(x, y, transpose(P), 
    xlabel = "x", ylabel = "y")
    #savefig("Reprezentari_Grafice\\Denumire.png")
end