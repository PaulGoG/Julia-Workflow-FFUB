#=
This part of the program is optional. It reads the output of the main program and processes it so that
it can be compared with experimental data (output data is averaged over experimental yield distributions)
=#

dY = CSV.read(yield_distribution_filename, DataFrame; delim = yield_distribution_delimiter, ignorerepeated = true, header = yield_distribution_header, skipto = yield_distribution_firstdataline)
println("reading $yield_distribution_filename done!")

#Outputs useful yield distribution formats for the program
function Process_yield_data(fragmdomain, dY)
    
end

if isassigned(dY.A[dY.A .< A₀/2], 1)
    println("A_L + A_H")
else
    println("A_H")
end