function DataDictionary(time_start,time_stop,time_step)

	# load stoichiometry
	stoichiometric_matrix = readdlm("./stoichiometry.dat")
	# Setup default flux bounds array
	default_bounds_array = [
		0 100.0; # 1 1.0*m_A_e<uptake:>1.0*m_A_c
		0 1.0; # 2 1.0*m_A_c<catalyze:>1.0*m_B_c
		0 1.0; # 3 1.0*m_B_c<catalyze:>1.0*m_A_c
		0 2.0; # 4 1.0*m_A_c<catalyze:>1.0*m_C_c
		0 100.0; # 5 1.0*m_C_c<catalyze:>1.0*m_B_c
		0 100.0; # 6 1.0*m_B_c<secrete:>1.0*m_B_e
		0 100.0; # 7 1.0*m_C_c<secrete:>1.0*m_C_e
	]

	# Setup default species bounds array
	species_bounds_array = [
		-10.0 10.0; # 1 m_A_e
		-10.0 10.0; # 2 m_B_e
		-10.0 10.0; # 3 m_C_e
		 0.0 0.0; # 4 m_A_c
		 0.0 0.0; # 5 m_B_c
		 0.0 0.0; # 6 m_C_c
	]

	# Setup the objective coefficient array
	objective_coefficient_array = [
		0.0; # 1 1.0*m_A_e<uptake:>1.0*m_A_c
		0.0; # 2 1.0*m_A_c<catalyze:>1.0*m_B_c
		0.0; # 3 1.0*m_B_c<catalyze:>1.0*m_A_c
		0.0; # 4 1.0*m_A_c<catalyze:>1.0*m_C_c
		0.0; # 5 1.0*m_C_c<catalyze:>1.0*m_B_c
		0.0; # 6 1.0*m_B_c<secrete:>1.0*m_B_e
		0.0; # 7 1.0*m_C_c<secrete:>1.0*m_C_e
	]

	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [
		"1.0*m_A_e<uptake:>1.0*m_A_c"
		"1.0*m_A_c<catalyze:>1.0*m_B_c"
		"1.0*m_B_c<catalyze:>1.0*m_A_c"
		"1.0*m_A_c<catalyze:>1.0*m_C_c"
		"1.0*m_C_c<catalyze:>1.0*m_B_c"
		"1.0*m_B_c<secrete:>1.0*m_B_e"
		"1.0*m_C_c<secrete:>1.0*m_C_e"
	]

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [
		"m_A_e"
		"m_B_e"
		"m_C_e"
		"m_A_c"
		"m_B_c"
		"m_C_c"
	]

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["is_minimum_flag"] = is_minimum_flag
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end



# maximize_product_dictionary
function maximize_product_dictionary(time_start,time_stop,time_step)
	# load the data dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)
	# Modify the data dictionary -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[6] = -1
	objective_coefficient_array[7] = -1
	# return -
	return data_dictionary
end
