#=
This part of the program is optional. It reads the output of the main program and processes it so that
it can be compared with experimental data (output data is averaged over experimental yield distributions)
=#

dY = CSV.read(yield_distribution_filename, DataFrame; delim = yield_distribution_delimiter, ignorerepeated = true, header = yield_distribution_header, skipto = yield_distribution_firstdataline)
println("reading $yield_distribution_filename done!")

#Outputs useful yield distribution formats for the program
function Process_yield_data(A_0, fragmdomain, dY)
    y_A_Z_TKE =  Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    if !isassigned(dY.A[dY.A .< A_0/2], 1)
        aux_dY = copy(dY)
        for A_H in unique(aux_dY.A)
            for TKE in unique(aux_dY.TKE[aux_dY.A .== A_H])
                if A_H != A_0 - A_H
                    push!(dY.A, A_0 - A_H)
                    push!(dY.TKE, TKE)
                    push!(dY.Value, dY.Value[(dY.A .== A_H) .& (dY.TKE .== TKE)][1])
                    push!(dY.σ, dY.σ[(dY.A .== A_H) .& (dY.TKE .== TKE)][1])
                end
            end
        end
    end

    #return y_A_Z_TKE
end

#Output all Necessary Y: Y(A,Z,TKE) for calculations and Y(A), Y(Z), Y(N), Y(TKE) for outputs and averaging