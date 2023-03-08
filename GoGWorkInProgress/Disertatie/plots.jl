#Function bodies and function calls for plotting data
if !isdir("plots/")
    mkdir("plots/")
end

using Plots
plotlyjs(size = (600, 800))
#####
function Grid_builder(x::AbstractVector, y::AbstractVector, z::AbstractVector)
    
end

#Testing for plots
Grid = zeros(length(unique(y_Ap_Z.A)), length(unique(y_Ap_Z.Z)))
i1 = 0
for A in unique(y_Ap_Z.A)
    i1 += 1
    i2 = 1
    for Z in unique(y_Ap_Z.Z)
        if isassigned(y_Ap_Z.Value[(y_Ap_Z.A .== A) .& (y_Ap_Z.Z .== Z)], 1)
            Grid[i1, i2] = y_Ap_Z.Value[(y_Ap_Z.A .== A) .& (y_Ap_Z.Z .== Z)][1]  
        end
        i2 += 1
    end
end

surface(Grid', zscale=:log10, zlims=(1e-1, maximum(Grid)))