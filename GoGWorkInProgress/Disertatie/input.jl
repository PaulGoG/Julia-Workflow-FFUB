file_output_identifier = "ID_TEST"                      #Identifier used in the filenames   
fissionant_nucleus_identifier = "CF52"                  #Symbol identifier of the fissionant nucleus
A₀, Z₀ = 252, 98                                        #(A,Z) of the fissionant nucleus
A_H_min, A_H_max = 126, 174                             #Heavy fission fragment range 
TKE_min, TKE_max, TKE_step = 130.0, 230.0, 2.0          #Total Kinetic Energy range and step in MeV
fission_type = "SF"                                     #SF or (n,f)
E_incident = 0.0                                        #Incident particle energy in MeV
density_parameter_type = "BSFG"                         #GC or BSFG
evaporation_cs_type = "VARIABLE"                        #CONSTANT or VARIABLE
isobaric_distribution_type = "DATA"                     #DATA or MEAN_VALUES
No_ZperA = 5                                            #Number of normally distributed Z considered for each A
txe_partitioning_type = "PARAM"                         #MSCZ, PARAM or RT
secondary_outputs = true                                #Yield-averaged outputs true or false selector