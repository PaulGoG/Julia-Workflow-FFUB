# Plotarea unei matrici NxN
function Reprezentari(x, y, Matrix)
    heatmap(x, y, transpose(Matrix), legend = false,
    xlabel = "x", ylabel = "y")
    contourf(x, y, transpose(Matrix), legend = false,
    xlabel = "x", ylabel = "y")
    surface(x, y, transpose(Matrix), legend = false,
    xlabel = "x", ylabel = "y")
    #savefig("Reprezentari_Grafice\\Denumire.png")
end