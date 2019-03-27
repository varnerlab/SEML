% FVA calling GLPK
% Solving the initial problem first, then does FVA
function [Flux_Array] = FVA(dataDictionary, MIN_MAX_FLAG)

	% Get some stuff from the DD -
	STM = dataDictionary('stoichiometric_matrix');
	[NUM_Species, NUM_Var] = size(STM);
	OBJVECTOR = dataDictionary('objective_coefficient_array');

	% Get Flux bounds from the DD -
	FluxBounds = dataDictionary('default_flux_bounds_array');
	FluxLB = FluxBounds(:,1);
	FluxUB = FluxBounds(:,2);
    % Formulate the VARTYPE vector (Flux)-
    VARTYPE = blanks(NUM_Var)';
    for id = 1:NUM_Var
        VARTYPE(id, 1) = 'C';
    end

	% Get species bounds
	SpeciesBounds = dataDictionary('species_bounds_array');
	% Get numbers
	NUM_Unbalanced = dataDictionary('extra_species_num');
	NUM_Balanced = NUM_Species - NUM_Unbalanced;

	% Equality constraints
	Aeq = STM((NUM_Unbalanced+1):NUM_Species, :);
	% Formulate the bV -
	bVEq = zeros(NUM_Balanced,1);

	% Inequality constraints
	UNBALANCED_STM = STM(1:NUM_Unbalanced, :);
	bVLB = SpeciesBounds(1:NUM_Unbalanced, 1);
	bVUB = SpeciesBounds(1:NUM_Unbalanced, 2);

    % Formulate the matrix of constraint coefficients
	A = [UNBALANCED_STM ; UNBALANCED_STM ; Aeq];
    % Formulate the column array containing the right-hand side value for
    % each constraint in the constraint matrix.
	bV = [bVUB ; bVLB; bVEq];
    % Formulate the CTYPE vector for each constraint
    CTYPE = blanks(length(bV))';
    for id = 1:NUM_Unbalanced
        CTYPE(id, 1) = 'U';
        CTYPE(id+NUM_Unbalanced, 1) = 'L';
    end
    offset = 2*NUM_Unbalanced;
    for id = (1+offset):(NUM_Balanced+offset)
        CTYPE(id) = 'S';
    end

    % ---------------------------------------------
    Flux_Array = zeros(NUM_Var, 3);

    % Phase One: Solve the initial problem P first -
	% Set the sense flag
    if (MIN_MAX_FLAG == -1)
		sense = -1;
    else
	 	sense = 1;
    end
	% Set some values for the optional parameters
	param.msglev = 3;
    param.tmlim = 3600;
	% Call the glpk solver -
	[FLOW,FVAL,P_status, P_extra] = glpk(OBJVECTOR, A, bV, FluxLB, FluxUB, CTYPE, VARTYPE, sense, param);
    if (P_status == 109)
        fprintf('failed to solve in 1h, try 10h...');
        param.tmlim = 3600*10;
        [FLOW,FVAL,P_status, P_extra] = glpk(OBJVECTOR, A, bV, FluxLB, FluxUB, CTYPE, VARTYPE, sense, param);
    end
    if (P_status == 109)
        fprintf('failed to find the optima after running for 10h');
    end
    for i = 1:NUM_Var
        Flux_Array(i, 1) = FLOW(i);
    end

    % Phase Two: FVA
    param.presol = 0;  % Turn off presolver
    new_OBJ = zeros(NUM_Var, 1);
    new_A = [A; OBJVECTOR];
    new_bV = [bV; FVAL];
    if (MIN_MAX_FLAG == -1)
        new_CT = [CTYPE; 'L'];
    else
        new_CT = [CTYPE; 'U'];
    end
    for round_id = 1:2
        if round_id == 1
            new_sense = 1;
        else
            new_sense = -1;
        end
        for i = 1:NUM_Var
            new_OBJ(i) = 1;
            param.tmlim = 3600;
            [F,fopt,S, E] = glpk(new_OBJ, new_A, new_bV, FluxLB, FluxUB, new_CT, VARTYPE, new_sense, param);
            if (S == 109)
                fprintf('failed to solve in 1h, try 10h...');
                param.tmlim = 3600*10;
                [F,fopt,S, E] = glpk(new_OBJ, new_A, new_bV, FluxLB, FluxUB, new_CT, VARTYPE, new_sense, param);
            end
            if (S == 109)
                fprintf('failed to find the optima after running for 10h');
            end
            Flux_Array(i, round_id + 1) = fopt;
            new_OBJ(i) = 0;
        end
    end

return;
