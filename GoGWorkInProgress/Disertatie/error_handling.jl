#=
Flow control and error handling for input parameters and filenames
=#
#####
include("input.jl")

if !isfile(mass_excess_filename)
    error("$mass_excess_filename does not exist or it cannot be accesed by the program!")
end

if fission_type != "SF" && fission_type != "(n,f)"
    error("$fission_type is not a valid input for fission_type!")
end

if density_parameter_type != "GC" && density_parameter_type != "BSFG"
    error("$density_parameter_type is not a valid input for density_parameter_type!")
end

density_parameter_filename