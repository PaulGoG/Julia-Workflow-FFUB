#Convert EXFOR x4 data format to a .dat file 
#dataindexoffset -> offset of first line of data below DATA tag in file
#selectedcols -> data columns that will be saved in the .dat file

cd(@__DIR__); cd("raw_data/")
using CSV, DataFrames

rawdatafile_name = "B.I.Starostov_4_1985"
dataindexoffset = 3
selectedcols = [1, 2, 3]

file_data = readlines(string(rawdatafile_name, ".txt"))
DATA_index = Int[]
for i in eachindex(file_data)
    if !isnothing(findfirst("DATA     ", file_data[i]))
        push!(DATA_index, i)
    end
end

DATA_index_start = first(DATA_index) + dataindexoffset
DATA_index_end = last(DATA_index) - 1
numerical_data = file_data[DATA_index_start:DATA_index_end]

aux_data_1 = filter.(x -> x != "", split.(numerical_data, " "))
aux_data_2 = [parse.(Float64, aux_data_1[i]) for i in eachindex(aux_data_1)]

if length(selectedcols) == length(aux_data_2[1])
    aux_data_3 = aux_data_2[1]
elseif length(selectedcols) > length(aux_data_2[1])
    for j in 1:(length(selectedcols) - length(aux_data_2[1]))
        push!(aux_data_2[1], 0.0)
    end
    aux_data_3 = aux_data_2[1]
else
    aux_data_3 = aux_data_2[1]
    deletecolsindex = deleteat!([i for i in eachindex(aux_data_3)], selectedcols)
    deleteat!(aux_data_3, deletecolsindex)
end

for i in 2:lastindex(aux_data_2)
    if length(selectedcols) == length(aux_data_2[i])
        aux_data_3 = hcat(aux_data_3, aux_data_2[i])
    elseif length(selectedcols) > length(aux_data_2[i])
        for j in 1:(length(selectedcols) - length(aux_data_2[i]))
            push!(aux_data_2[i], 0.0)
        end
        aux_data_3 = hcat(aux_data_3, aux_data_2[i])
    else 
        deletecolsindex = deleteat!([j for j in eachindex(aux_data_2[i])], selectedcols)
        deleteat!(aux_data_2[i], deletecolsindex)
        aux_data_3 = hcat(aux_data_3, aux_data_2[i])
    end
end
aux_data_3 = aux_data_3'

if length(selectedcols) > length(aux_data_3[1, :])
    aux_data_4 = DataFrame(aux_data_3, :auto)[:, selectedcols]
else 
    aux_data_4 = DataFrame(aux_data_3, :auto)
    selectedcols = [i for i in eachindex(selectedcols)]
end


aux_data = DataFrame(
                    E = aux_data_4[!, selectedcols[1]], 
                    ratioMxwll = aux_data_4[!, selectedcols[2]],
                    erro = A = aux_data_4[!, selectedcols[3]]
                    )
display(aux_data)

cd(@__DIR__); cd("output/")
CSV.write(string(rawdatafile_name, ".dat"), aux_data, delim=' ')