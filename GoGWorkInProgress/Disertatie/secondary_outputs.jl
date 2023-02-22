#=
This program reads the output of the main program and processes it so that
it can be compared with experimental data and used for primary validations (raw output data)
or secondary validations (output data averaged over experimental yield distributions)
=#

include("error_handling.jl")

#Yield distribution filename
yield_distribution_filename = "U5YATKE.SRE"


if !isfile(yield_distribution_filename)
    error("$yield_distribution_filename does not exist at input_data/  PATH!")
end