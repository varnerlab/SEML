from __future__ import division
from FluxDriver import FluxDriver
from DataDictionary import maximize_product_dictionary
from FVA import calculate_flux_variabilty

# load the data dictionary -
data_dictionary = maximize_product_dictionary(0,10,1)

# solve the lp problem -
objective_value, flux_array, uptake_array, exit_flag = FluxDriver(data_dictionary)
calculated_flux_array = calculate_flux_variabilty(data_dictionary)

print "objective_value: ", objective_value
print "flux array: ", flux_array
print "uptake array: ", uptake_array
print "FVA results: \n", calculated_flux_array
