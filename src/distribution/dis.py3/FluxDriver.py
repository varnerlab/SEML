from scipy.optimize import linprog
import numpy as np
import sys

def FluxDriver(dataDictionary):
    # Get some stuff from the DF -
    STM = dataDictionary["stoichiometric_matrix"]
    OBJVECTOR = dataDictionary["objective_coefficient_array"]
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

    # Call the LP solver -
    opt = {"disp": True, "tol": 1.0E-8}
    res = linprog(OBJVECTOR, A_ub=A, b_ub=bV, A_eq=Aeq, b_eq=bVEq,
                  bounds=FluxBounds, options=opt)
    if res.status == 1:
        opt = {"disp": True, "tol": 1.0E-8, "maxiter": 10000}
        res = linprog(OBJVECTOR, A_ub=A, b_ub=bV, A_eq=Aeq, b_eq=bVEq,
                      bounds=FluxBounds, options=opt)
    if res.status == 2:
        sys.exit()

    UPTAKE = np.dot(UNBALANCED_STM, res.x)

    return res.fun, res.x, UPTAKE, res.status
