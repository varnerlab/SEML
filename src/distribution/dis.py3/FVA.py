from scipy.optimize import linprog
import numpy as np
import sys

# FVA in python
def calculate_flux_variabilty(dataDictionary):
    # Get some stuff from the DF -
    STM = dataDictionary["stoichiometric_matrix"]
    OBJVECTOR = dataDictionary["objective_coefficient_array"]
    number_of_fluxes = np.size(OBJVECTOR)
    # Get Flux bounds sequence from the DF -
    FluxBounds = dataDictionary["default_flux_bounds_array"]
    # Get species bounds
    SpeciesBounds = dataDictionary["species_bounds_array"]
    # Get numbers
    NUM_Unbalanced = dataDictionary["extra_species_num"]
    NUM_Speices = len(SpeciesBounds)
    NUM_Balanced = NUM_Speices - NUM_Unbalanced

    # Equality constraints
    Aeq = STM[NUM_Unbalanced:NUM_Speices]
    # Formulate the bV -
    bVEq = np.zeros(NUM_Balanced)

    # Inequality constraints
    UNBALANCED_STM = STM[0:NUM_Unbalanced]
    bVLB = [-SpeciesBounds[i][0] for i in range(NUM_Unbalanced)]
    bVUB = [SpeciesBounds[i][1] for i in range(NUM_Unbalanced)]
    A = np.concatenate((UNBALANCED_STM, -1*UNBALANCED_STM))
    bV = np.concatenate((bVUB, bVLB))

    # initialize the output array -
    calculated_flux_array = np.zeros((np.size(OBJVECTOR), 3))

    # Phase 1: solve the original problem as a minimization problem-
    # default parameter setting of linprog: maxiter=1000, disp=False, tol=1.0E-12,
    opt = {"disp": True, "tol": 1.0E-8}
    FBA_result = linprog(OBJVECTOR, A_ub=A, b_ub=bV, A_eq=Aeq, b_eq=bVEq,
                 bounds=FluxBounds, options=opt)
    if FBA_result.status == 1:
        print("\nrerun linprog with 10*maxiter")
        opt = {"disp": True, "tol": 1.0E-8, "maxiter": 10000}
        FBA_result = linprog(OBJVECTOR, A_ub=A, b_ub=bV, A_eq=Aeq, b_eq=bVEq,
                     bounds=FluxBounds, options=opt)
    if FBA_result.status == 2:
        sys.exit()
    for i in range(number_of_fluxes):
        calculated_flux_array[i, 0] = FBA_result.x[i]

    # Phase 2: solve the FVA problem
    new_A = np.concatenate((A, [OBJVECTOR]))
    new_bV = np.concatenate((bV, [FBA_result.fun]))
    new_OBJVECTOR = np.zeros(number_of_fluxes)
    # main FVA loop -
    for round_index in range(2):
        # min or max
        if round_index == 0:
            OptDir = 1
        else:
            OptDir = -1
        for i in range(number_of_fluxes):
            new_OBJVECTOR[i] = OptDir
            opt = {"disp": True, "tol": 1.0E-8}
            res = linprog(new_OBJVECTOR, A_ub=new_A, b_ub=new_bV, A_eq=Aeq,
                          b_eq=bVEq, bounds=FluxBounds, options=opt)
            if res.status == 1:
                opt = {"disp": True, "tol": 1.0E-8, "maxiter": 10000}
                res = linprog(new_OBJVECTOR, A_ub=new_A, b_ub=new_bV, A_eq=Aeq,
                              b_eq=bVEq, bounds=FluxBounds, options=opt)
            if res.status != 2:
                calculated_flux_array[i, round_index+1] = OptDir*res.fun
            else:
                print("ERROR: Unable to solve modified problem for flux indexed as", i)
            new_OBJVECTOR[i] = 0

    return calculated_flux_array
