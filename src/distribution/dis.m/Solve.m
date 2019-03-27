clc;
% load the data dictionary -
data_dictionary = maximizeProductDictionary(0,10,1);

% solve the lp problem -
[Flux,fVal,UptakeRate,EXITFLAG] = FluxDriver(data_dictionary, 1);
calculated_flux_array = FVA(data_dictionary, 1);

fprintf('objective value is %d \n', fVal);
fprintf('flux array: \n');
disp(Flux);
fprintf('uptake rate: \n');
disp(UptakeRate);
fprintf('FVA results: \n');
disp(calculated_flux_array);
