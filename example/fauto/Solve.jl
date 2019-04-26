include("Include.jl")

# load the data dictionary -
data_dictionary = maximize_product_dictionary(0,10,1)

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
(calculated_flux_array, dual_value_array) = calculate_flux_variabilty(data_dictionary)

# FBA
println("objective: $objective_value")
println("flux array: ", flux_array)
println("uptake array: ", uptake_array)

# FVA
println("FVA results: ")
for i = 1:size(calculated_flux_array)[1]
  println("    ", calculated_flux_array[i,:])
end
