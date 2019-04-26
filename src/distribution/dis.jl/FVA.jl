function calculate_flux_variabilty(data_dictionary)

    # --- Phase 0: GET SOME STUFF ---------------------------------------------------#
    # Get the stoichiometric_matrix from data_dictionary -
    stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]
    (number_of_species,number_of_fluxes) = size(stoichiometric_matrix);

    # initialize the output array -
    calculated_flux_array = zeros(number_of_fluxes,3)

    # set some constants -
    TIME_RESTART_LIM = 60
    EXIT_FLAG_SUCCESS = 0
    TOL = 1e-6
    # ------------------------------------------------------------------------------- #

    # --- Phase 1: SOLVE THE ORIGINAL PROBLEM --------------------------------------- #

    # Setup the GLPK problem -
    lp_problem = GLPK.Prob();
    GLPK.set_prob_name(lp_problem, "sample");
    GLPK.set_obj_name(lp_problem, "objective")

    # min -or- max?
    min_flag = data_dictionary["is_minimum_flag"];
    if min_flag == true
    	GLPK.set_obj_dir(lp_problem, GLPK.MIN);
    else
    	GLPK.set_obj_dir(lp_problem, GLPK.MAX);
    end

    # Set the number of constraints and fluxes -
    GLPK.add_rows(lp_problem, number_of_species);
    GLPK.add_cols(lp_problem, number_of_fluxes);

    # Set solver parameters
    solver_parameters = GLPK.SimplexParam();
    solver_parameters.msg_lev = GLPK.MSG_ERR;
    solver_parameters.presolve = GLPK.ON;

    # Setup flux bounds -
    default_bounds_array = data_dictionary["default_flux_bounds_array"]
    (number_of_fluxes,number_of_bounds) = size(default_bounds_array)
    for flux_index = 1:number_of_fluxes

    	flux_lower_bound = default_bounds_array[flux_index,1]
    	flux_upper_bound = default_bounds_array[flux_index,2]

    	# Check bounds type ... default is DB -
    	if (flux_upper_bound == flux_lower_bound)
    		flux_constraint_type = GLPK.FX
    	else
    		flux_constraint_type = GLPK.DB
    	end

    	# flux symbol? (later use name - for now, fake it)
    	flux_symbol = "R_"*string(flux_index)

    	# Set the bounds in GLPK -
    	GLPK.set_col_name(lp_problem, flux_index, flux_symbol);
    	GLPK.set_col_bnds(lp_problem, flux_index, flux_constraint_type, flux_lower_bound, flux_upper_bound);
    end

    # Setup objective function -
    objective_coefficient_array = data_dictionary["objective_coefficient_array"]
    for (flux_index,obj_coeff) in enumerate(objective_coefficient_array)
    	# Set the objective function value in GLPK -
    	GLPK.set_obj_coef(lp_problem, flux_index, obj_coeff);
    end

    # Setup problem constraints for the metabolites -
    species_bounds_array = data_dictionary["species_bounds_array"]
    for species_index = 1:number_of_species
    	species_lower_bound = species_bounds_array[species_index,1]
    	species_upper_bound = species_bounds_array[species_index,2]
    	# defualt
    	species_constraint_type = GLPK.FX
    	if (species_lower_bound != species_upper_bound)
    		species_constraint_type = GLPK.DB
    	end
    	# set the symbol -
    	species_symbol = "x_"*string(species_index)
    	# Set the species bounds in GLPK -
    	GLPK.set_row_name(lp_problem, species_index, species_symbol);
    	GLPK.set_row_bnds(lp_problem, species_index, species_constraint_type, species_lower_bound, species_upper_bound);
    end

    # Setup the stoichiometric array -
    counter = 1;
    row_index_array = zeros(Int,number_of_species*number_of_fluxes);
    col_index_array = zeros(Int,number_of_species*number_of_fluxes);
    species_index_vector = collect(1:number_of_species);
    flux_index_vector = collect(1:number_of_fluxes);
    flat_stoichiometric_array = zeros(Float64,number_of_species*number_of_fluxes);
    for species_index in species_index_vector
    	for flux_index in flux_index_vector
    		row_index_array[counter] = species_index;
    		col_index_array[counter] = flux_index;
    		flat_stoichiometric_array[counter] = stoichiometric_matrix[species_index,flux_index];
    		counter+=1;
    	end
    end

    # Passed in the matrix
    GLPK.load_matrix(lp_problem, number_of_species*number_of_fluxes, row_index_array, col_index_array, flat_stoichiometric_array);

    # Call the solver -
    exit_flag = GLPK.simplex(lp_problem, solver_parameters);

    # Get the objective function value -
    objective_value = GLPK.get_obj_val(lp_problem);

    # Get the calculated flux values for the optimal distrubution from GLPK -
    for flux_index in flux_index_vector
    	calculated_flux_array[flux_index,1] = GLPK.get_col_prim(lp_problem, flux_index);
    end

    # Get the dual values -
    dual_value_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
    	dual_value_array[flux_index] = GLPK.get_col_dual(lp_problem, flux_index);
    end
    # ------------------------------------------------------------------------------- #

    # --- Phase 2: SOLVE THE FVA PROBLEM -------------------------------------------- #
    # retrieve the optimal value
    target_value = 0.0
    if (GLPK.get_obj_dir(lp_problem) == GLPK.MAX)
        target_value = floor(objective_value/TOL)*(TOL);
    else
        target_value = ceil(objective_value/TOL)*(TOL);
    end

    # add the constraint from phase I (the object value) -
    m = GLPK.add_rows(lp_problem,1)
    if (GLPK.get_obj_dir(lp_problem) == GLPK.MAX)
        GLPK.set_row_bnds(lp_problem, m, GLPK.LO,target_value,0.0);
    else
        GLPK.set_row_bnds(lp_problem, m, GLPK.UP,0.0,target_value);
    end
    index_array = Int[]
    value_array = Float64[]
    for flux_index = 1:number_of_fluxes
        push!(index_array,flux_index)
        push!(value_array,GLPK.get_obj_coef(lp_problem, flux_index))
    end
    GLPK.set_mat_row(lp_problem,m,number_of_fluxes,index_array,value_array)

    # zero out the objective coefficients -
    for flux_index = 1:number_of_fluxes
        GLPK.set_obj_coef(lp_problem, flux_index, 0.0);
    end

    # set the parameters
    solver_parameters.presolve = GLPK.OFF;
    solver_parameters.msg_lev = GLPK.MSG_OFF;
    solver_parameters.tm_lim = 1000*TIME_RESTART_LIM;

    # main FVA loop -
    for round_index = 1:2
        # min or maximizing the flux?
        if round_index == 1
            GLPK.set_obj_dir(lp_problem, GLPK.MIN);
        else
            GLPK.set_obj_dir(lp_problem, GLPK.MAX);
        end
        # run for each flux
        for flux_index = 1:number_of_fluxes
            # run the calc -
            GLPK.set_obj_coef(lp_problem, flux_index, 1.0);
            # ??? so by default, it would be fastFVA?
            exit_flag = GLPK.simplex(lp_problem, solver_parameters);
            # try to retrieve results
            if (exit_flag != EXIT_FLAG_SUCCESS)
                solver_parameters.tm_lim = 10000*TIME_RESTART_LIM;
                solver_parameters.presolve = GLPK.ON
                GLPK.adv_basis(lp_problem, 0);  # not fastFVA
                exit_flag = GLPK.simplex(lp_problem, solver_parameters);
                if (exit_flag != EXIT_FLAG_SUCCESS)
                    # ok, something muy malo is happening. throw an error ...
                    error("ERROR: Unable to solve modified problem. Exiting ... ")
                end
                # reset the problem options -
                solver_parameters.presolve = GLPK.OFF;
                solver_parameters.tm_lim = 1000*TIME_RESTART_LIM;
            end
            # reset the objective coefficient back to zero -
            GLPK.set_obj_coef(lp_problem,flux_index,0.0)
            # store the results
            if (GLPK.get_obj_dir(lp_problem) == GLPK.MIN)
                calculated_flux_array[flux_index,2] = GLPK.get_obj_val(lp_problem)
            else
                calculated_flux_array[flux_index,3] = GLPK.get_obj_val(lp_problem)
            end
        end
    end
    # ------------------------------------------------------------------------------- #

    # return the array -
    return (calculated_flux_array,dual_value_array)
end
