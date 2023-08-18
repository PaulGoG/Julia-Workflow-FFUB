#Additional options for reading external files and true/false selector variables 
########################################################################################################

mass_excess_filename = "AME2020.ANA"
mass_excess_header = ["Z", "A", "Symbol", "D", "σ_D"]
mass_excess_delimiter = ' '
mass_excess_firstdataline = 1

density_parameter_filename = "SZSN.GC"
density_parameter_header = ["n", "S_Z", "S_N"]
density_parameter_delimiter = ' '
density_parameter_firstdataline = 2

isobaric_distribution_filename = "DeltaZ_rms_A.$fissionant_nucleus_identifier"
isobaric_distribution_header = ["A", "ΔZ_A", "rms_A"]
isobaric_distribution_delimiter = ' '
isobaric_distribution_firstdataline = 2

txe_partitioning_filename = ""
txe_partitioning_header = ["A", "Z", "Value"]
txe_partitioning_delimiter = ' '
txe_partitioning_firstdataline = 2

txe_partitioning_segmentpoints = [(126, 0.5), (130, 0.14), (150, 0.607), (174, 0.865)]

yield_distribution_filename = "$(fissionant_nucleus_identifier)YATKE.VES"
yield_distribution_header = ["A", "TKE", "Value", "σ"]
yield_distribution_delimiter = ' '
yield_distribution_firstdataline = 2

secondary_output_Yield = false
secondary_output_TXE_Q = false
secondary_output_E_excitation = false
secondary_output_nu = false
secondary_output_Ap = false
secondary_output_T = false
secondary_output_Tₖ = false
secondary_output_avg_ε = false
secondary_output_avg_εₖ = false
secondary_output_Eᵣ = false

ΔT, ΔTₖ, Δavg_ε, Δavg_εₖ, ΔEᵣ = 1e-2, 1e-2, 5e-2, 5e-2, 5e-1

neutron_spectrum = false
E_min, E_max, E_step = 1e-6, 20.0, 5e-2

generate_plots = false
multidim_plots = false
resolution_scale = 250.0
aspect_ratio = (3, 2)