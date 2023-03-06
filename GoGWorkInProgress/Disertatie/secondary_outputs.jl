#=
This part of the program is optional. It reads the output of the main program and processes it so that
it can be compared with experimental data (output data is averaged over experimental yield distributions)
=#

dY = CSV.read(yield_distribution_filename, DataFrame; delim = yield_distribution_delimiter, ignorerepeated = true, header = yield_distribution_header, skipto = yield_distribution_firstdataline)
println("reading $yield_distribution_filename done!")

#Compute Y(A,Z,TKE) form experimental Y(A,TKE) data
function Process_yield_data(A_0, fragmdomain, dY)
    y_A_Z_TKE =  Distribution(Int[], Int[], Float64[], Int[], Float64[], Float64[])
    #Completes fragmentation domain in case data provided only for HF
    if !isassigned(dY.A[dY.A .< A_0/2], 1)
        for A in unique(sort(dY.A, rev=true))
            for Z in fragmdomain.Z[fragmdomain.A .== A_0 - A]
                P_A_Z = fragmdomain.Value[(fragmdomain.A .== A_0 - A) .& (fragmdomain.Z .== Z)][1]
                for TKE in unique(sort(dY.TKE[dY.A .== A]))
                    val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                    if A != A_0 - A
                        push!(y_A_Z_TKE.A, A_0 - A)
                        push!(y_A_Z_TKE.Z, Z)
                        push!(y_A_Z_TKE.TKE, TKE)
                        push!(y_A_Z_TKE.Value, val)
                        push!(y_A_Z_TKE.σ, σ)
                    end
                end
            end
        end
    end
    for A in unique(sort(dY.A))
        for Z in fragmdomain.Z[fragmdomain.A .== A]
            P_A_Z = fragmdomain.Value[(fragmdomain.A .== A) .& (fragmdomain.Z .== Z)][1]
            for TKE in unique(sort(dY.TKE[dY.A .== A]))
                val = P_A_Z * dY.Value[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                σ = P_A_Z * dY.σ[(dY.A .== A) .& (dY.TKE .== TKE)][1]
                if A != A_0 - A
                    push!(y_A_Z_TKE.A, A)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                elseif !isassigned(y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)], 1)
                    push!(y_A_Z_TKE.A, A)
                    push!(y_A_Z_TKE.Z, Z)
                    push!(y_A_Z_TKE.TKE, TKE)
                    push!(y_A_Z_TKE.Value, val)
                    push!(y_A_Z_TKE.σ, σ)
                else 
                    y_A_Z_TKE.Value[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)] .+= val
                    y_A_Z_TKE.σ[(y_A_Z_TKE.A .== A) .& (y_A_Z_TKE.Z .== Z) .& (y_A_Z_TKE.TKE .== TKE)] .= σ*sqrt(2)
                end
            end
        end
    end
    #Renormalizes data
    f = 200/sum(y_A_Z_TKE.Value)
    y_A_Z_TKE.Value .= y_A_Z_TKE.Value .* f
    y_A_Z_TKE.σ .= y_A_Z_TKE.σ .* f
    return y_A_Z_TKE
end
#Average f(A,Z,TKE) over p(Z,A)
#Average f(A,Z,TKE) over Y(A,Z,TKE) so it becomes f(A)
#Average f(A,Z,TKE) over Y(A,Z,TKE) so it becomes f(TKE)
#Average f(A,Z,TKE) over Y(A,Z,TKE) to get average value <f>
#Obtain Y(Aₚ) distribution from Y(A) & ν(A)
#Obtain Y(Z, Aₚ) distribution from Y(A,Z) & ν(A,Z)