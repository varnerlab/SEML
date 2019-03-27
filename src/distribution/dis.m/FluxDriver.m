% FluxDriver.m calling GLPK
function [FLOW,FVAL,UPTAKE,status] = FluxDriver(dataDictionary, MIN_MAX_FLAG)

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
	[FLOW,FVAL,status, extra] = glpk(OBJVECTOR, A, bV, FluxLB, FluxUB, CTYPE, VARTYPE, sense, param);
    if (status == 109)
        fprintf('failed to solve in 1h, try 10h...');
        param.tmlim = 3600*10; 
        [FLOW,FVAL,status, extra] = glpk(OBJVECTOR, A, bV, FluxLB, FluxUB, CTYPE, VARTYPE, sense, param);
    end 
    if (status == 109) 
        fprintf('failed to find the optima after running for 10h');
    end 

	UPTAKE = UNBALANCED_STM*FLOW;
return;
