# Plotarea unei matrici NxN
function Reprezentari(x, y, P)
    heatmap(x, y, transpose(P),  
    xlabel = "x", ylabel = "y")
    contourf(x, y, transpose(P), 
    xlabel = "x", ylabel = "y")
    surface(x, y, transpose(P),  
    xlabel = "x", ylabel = "y")
    #savefig("Reprezentari_Grafice\\Denumire.png")
end