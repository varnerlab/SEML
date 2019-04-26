# MATLAB

# specific to Julia, that is, related to the output file format (language).

#=
F: read the file containing constants
I:
O:
BA
=#
function m_include_constants_from_literature(src_file_name, pad_string)
  # create src_buffer -
  src_buffer = ""
  # path to distrubtion -
  path_to_src_file = src_file_name
  open(path_to_src_file,"r") do src_file
    for line in eachline(src_file)
      src_buffer *= pad_string*"$line"
    end
  end
  return src_buffer
end

#=
F: generate kinetic buffer
I:
O: return monod affinity constant symbol, mRNA species array,
  protein species array, W_string_array, disassociation_const_string_array.
BA
=#
function m_build_kinetics_buffer(all_species_dict::Dict{String, Int},
  all_rnx_list::Array, all_txtl_dict::Dict, sys2user::Dict)
  kinetics = "function [rnx_rate_vector, TX_rate_vector, TL_rate_vector] = Kinetics(X, data_dictionary)"
            #  "\n\n\t% load all species dictionary" *
            #  "\n\tall_species_dict = data_dictionary(\'all_species_dict\');  % String-->Int"
  # for signaling
  kinetics *= "\n\n\t% generate kinetics for signaling, assume saturation kinetics"
  totalRnxNo = length(all_rnx_list)
  kinetics *= "\n\t% load Monod affinity constants" *
              "\n\tMonodK = data_dictionary(\'MonodAffinityConstantDict\');  % String-->Float"
  MonodAffinityConstant_String_Array = Array{String,1}()  # collect monod affinity constant symbol
  kinetics *= "\n\t% load reaction kinetic constants" *
              "\n\tkcat = data_dictionary(\'kcat_signaling\');  % in order of rnx" *
              "\n\n\t% write reaction rate equations" *
              "\n\trnx_rate_vector = zeros($totalRnxNo, 1);"
  for (index, rnx) in enumerate(all_rnx_list)  # go thru every rnx
    tmp_rnx_rate = "kcat($index)"
    if isdefined(rnx, :catalysts)  # enzyme
      for token in rnx.catalysts
        tmp_rnx_rate *= "*X($(all_species_dict[token.oriBioName]))"
      end
    end
    if isdefined(rnx,:reactants)  # substrate
      for token in rnx.reactants
        tmp_rnx_rate *= "*X($(all_species_dict[token.oriBioName]))/(MonodK(\'MonodK~rnx$(index)~$(token.oriBioName)\')+X($(all_species_dict[token.oriBioName])))"
        push!(MonodAffinityConstant_String_Array, "MonodK~rnx$(index)~$(token.oriBioName)")
      end
    end
    kinetics *= "\n\trnx_rate_vector($index) = "*tmp_rnx_rate * ";"
  end


  # for TXTL
  W_string_array = Array{String,1}()  # W: weight of protein action on transcription
  disassociation_const_string_array = Array{String,1}()  # disassociation constant in tranfer function
  kinetics *= "\n\n\t% generate kinetics & control terms for TXTL;" *
              "\n\t% load data from data dictionary;" *
              "\n\tW_value_dict = data_dictionary(\'W_value_dict\');  % String-->Float" *
              "\n\tcoop = data_dictionary(\'cooperativity\');" *
              "\n\tdisassociation_constant_dict = data_dictionary(\'transferFunctionDisassociationConstantDict\');" *
              "\n\tbackground_control_term = data_dictionary(\'backgroundGeneExpressionControlTermDict\');" *
              "\n\tkcatTX = data_dictionary(\'kcatTranscription\');" *
              "\n\tRNAPconc = data_dictionary(\'RNAPConcentration\');" *
              "\n\tgeneConc = data_dictionary(\'avgGeneConcentration\');" *
              "\n\tsaturationTX = data_dictionary(\'transcriptionSaturationConstant\');" *
              "\n\tkcatTL = data_dictionary(\'kcatTranslation\');" *
              "\n\tRIBOconc = data_dictionary(\'RIBOConcentration\');" *
              "\n\tsaturationTL = data_dictionary(\'translationSaturationConstant\');"

  totalTXTLNo = length(all_txtl_dict)
  kinetics *= "\n\n\t% write control terms, transcription rate equations, and translation rate equations" *
             "\n\tTX_control_term = containers.Map();  % control term" *
             "\n\tTX_rate_vector = containers.Map();  % transcription rate" *
             "\n\tTL_rate_vector = containers.Map();  % translation rate"
  for (key, txtl) in all_txtl_dict # go thru every txtl
    tmp_protein_string = replace(key, sys2user["MRNA"] => sys2user["PROTEIN"], count=1)
    kinetics *= "\n\t% $key and $tmp_protein_string"
    up_factors_array = Array{String,1}()  # collection of upregulation factors name: targetedmRNA_factor(s)
    if !isempty(txtl.activationProtein)  # upregulation
      for token_array in txtl.activationProtein
        tmp_up = ""
        tmp_up_name = key*"_"
        for token in token_array
          tmp_W = "W~$(token.oriBioName)~$(key)"
          tmp_up *= "W_value_dict(\'$(tmp_W)\') * X($(all_species_dict[token.oriBioName]))^coop/((disassociation_constant_dict(\'KD~$(token.oriBioName)\'))^coop+X($(all_species_dict[token.oriBioName]))^coop)*"
          push!(W_string_array, tmp_W)
          push!(disassociation_const_string_array, "KD~$(token.oriBioName)")
          tmp_up_name *= token.oriBioName*"_"
        end
        tmp_up = chop(tmp_up)
        tmp_up_name = chop(tmp_up_name)
        push!(up_factors_array, tmp_up_name)
        kinetics *= "\n\t$(tmp_up_name) = $(tmp_up);"
      end
    end
    down_factors_array = Array{String,1}()  # collection of downregulation factors name: targetedmRNA_factor(s)
    if !isempty(txtl.inhibitionProtein)  # downregulation
      for token_array in txtl.inhibitionProtein
        tmp_down = ""
        tmp_down_name = key*"_"
        for token in token_array
          tmp_W = "W~$(token.oriBioName)~$(key)"
          tmp_down *= "W_value_dict(\'$(tmp_W)\') * X($(all_species_dict[token.oriBioName]))^coop/((disassociation_constant_dict(\'KD~$(token.oriBioName)\'))^coop+X($(all_species_dict[token.oriBioName]))^coop)*"
          push!(W_string_array, tmp_W)
          push!(disassociation_const_string_array, "KD~$(token.oriBioName)")
          tmp_down_name *= token.oriBioName*"_"
        end
        tmp_down = chop(tmp_down)
        tmp_down_name = chop(tmp_down_name)
        push!(down_factors_array, tmp_down_name)
        kinetics *= "\n\t$(tmp_down_name) = $(tmp_down);"
      end
    end
    # combine to form up and down, then merge up and down to form control term
    if !isempty(up_factors_array) && !isempty(down_factors_array)
      up_action = ""
      for up in up_factors_array
        up_action *= up*" +"
      end
      kinetics *= "\n\tup_action_on_$key = $(chop(up_action));"
      down_action = ""
      for down in down_factors_array
        down_action *= down*" +"
      end
      kinetics *= "\n\tdown_action_on_$key = $(chop(down_action));" *
                  "\n\tTX_control_term(\'$key\') = (background_control_term(\'$key\') + up_action_on_$key)/(1 + background_control_term(\'$key\') + up_action_on_$key + down_action_on_$key);"
    elseif !isempty(up_factors_array)
      up_action = ""
      for up in up_factors_array
        up_action *= up*" +"
      end
      kinetics *= "\n\tup_action_on_$key = $(chop(up_action));" *
                  "\n\tTX_control_term(\'$key\') = (background_control_term(\'$key\') + up_action_on_$key)/(1 + background_control_term(\'$key\') + up_action_on_$key);"
    elseif !isempty(down_factors_array)
      down_action = ""
      for down in down_factors_array
        down_action *= down*" +"
      end
      kinetics *= "\n\tdown_action_on_$key = $(chop(down_action));" *
                  "\n\tTX_control_term(\'$key\') = background_control_term(\'$key\')/(1 + background_control_term(\'$key\') + down_action_on_$key);"
    else
      # seems impossible, a TXTL dict is constructed b/c something acts on the key
    end
    kinetics *= "\n\tTX_rate_vector(\'$key\') = TX_control_term(\'$key\') * kcatTX * RNAPconc * geneConc / (saturationTX + geneConc);" *
                "\n\tTL_rate_vector(\'$tmp_protein_string\') = kcatTL * RIBOconc * X($(all_species_dict[key])) / (saturationTL + X($(all_species_dict[key])));"
  end

  kinetics *= "\n\nend"
  return (kinetics, MonodAffinityConstant_String_Array, W_string_array, disassociation_const_string_array)
end

function m_build_data_dictionary_buffer(host_type::AbstractString, all_species_array::Array,
  all_species2index_dict::Dict,
  rnx_species_array::Array, all_rnx_list::Array, all_txtl_dict::Dict,
  Monod_affinity_constant_array::Array, W_string_array::Array,
  disassociation_const_array::Array, mRNA_species_array::Array, protein_species_array::Array)
  # data dictionary
  buffer = "function dataDictionary = DataDictionary()"

  # stoichiometry
  buffer *= "\n\n\t% load stoichiometry"
  buffer *= "\n\tstoichiometric_matrix = load(\'stoichiometry.dat\');"
  No_species = length(all_species_array)
  # initial concentration
  buffer *= "\n\n\t% initial species concentration" *
    "\n\tinitial_condition = [\n\t\t"
  buffer *= join(["0.001, ...  % $i $(all_species_array[i])" for i = 1:(No_species - 1)], "\n\t\t")
  buffer *= "\n\t\t0.001  ...  % $No_species $(all_species_array[No_species])"
  buffer *= "\n\t];"
  # all species dictionary
  buffer *= "\n\n\t% all species dictionary"
  buffer *= "\n\tall_species_dict = containers.Map();  % String --> Int"
  #all_species_dict = Dict{String, Int64}()  # create all species dict: string --> int64
  for (st, index) in all_species2index_dict
    buffer *= "\n\tall_species_dict(\'$st\') = $index;"
  end
  # all species dictionary_reversed
  buffer *= "\n\n\t% all species dictionary_reversed"
  buffer *= "\n\tkey = keys(all_species_dict);" *
            "\n\tval = values(all_species_dict);" *
            "\n\tall_species_reversed_dict = containers.Map(val, key);  % Int --> String"

  # all rnx species array
  buffer *= "\n\n\t% all rnx species array"
  buffer *= "\n\trnx_species_array = {..." *
    join(["\n\t\t\'$st\'" for st in rnx_species_array], ", ...") *
    "  ...\n\t};"
  # all mRNA species array
  buffer *= "\n\n\t% all mRNA species array"
  buffer *= "\n\tmRNA_species_array = {..." *
    join(["\n\t\t\'$st\'" for st in mRNA_species_array], ", ...") *
    "  ...\n\t};"
  # all protein species array
  buffer *= "\n\n\t% all protein species array"
  buffer *= "\n\tprotein_species_array = {..." *
    join(["\n\t\t\'$st\'" for st in protein_species_array], ", ...") *
    "  ...\n\t};"
  # signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)
  buffer *= "\n\n\t% signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)"
  buffer *= "\n\tkcat_signaling = 1e-3*ones(1, $(length(all_rnx_list)));  % kcat[#reaction]: reaction name"
  for (id, rnx) in enumerate(all_rnx_list)
    buffer*= "\n\tkcat_signaling($id) = 1.0;  % kcat: $(rnx.rnxName)"
  end
  # Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate
  buffer *= "\n\n\t% Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate"
  buffer *= "\n\tMonodAffinityConstantDict = containers.Map();"
  for MK in Monod_affinity_constant_array
    # buffer *= "\n\t$MK = 1.0"
    # buffer *= "\n\tMonodAffinityConstant_dict(\'$MK\') = $MK"
    buffer *= "\n\tMonodAffinityConstantDict(\'$MK\') = 1.0;"
  end
  # W_value in transcription control term, nomenclature: W_targetmRNA_actor
  buffer *= "\n\n\n\t% W_value in transcription control term, nomenclature: W_targetmRNA_actor"
  buffer *= "\n\tW_value_dict = containers.Map();"
  for WS in W_string_array
    buffer *= "\n\tW_value_dict(\'$WS\') = 1.0;"
  end
  # disassociation constants in transfer function, nomenclature: KD_speciesName
  buffer *= "\n\n\t% disassociation constants in transfer function, nomenclature: KD_speciesName"
  buffer *= "\n\ttransferFunctionDisassociationConstantDict = containers.Map();"
  for DC in disassociation_const_array
    buffer *= "\n\ttransferFunctionDisassociationConstantDict(\'$DC\') = 1.0;"
  end
  # transcription specific correction factor
  buffer *= "\n\n\t% transcription specific correction factor"
  buffer *= "\n\ttranscriptionSpecificCorrectionFactorDict = containers.Map();"
  for mrna in mRNA_species_array
    buffer *= "\n\ttranscriptionSpecificCorrectionFactorDict(\'$mrna\') = 1.0;  % $mrna"
  end
  # background gene expression control term
  buffer *= "\n\n\t% background gene expression control term"
  buffer *= "\n\tbackgroundGeneExpressionControlTermDict = containers.Map();"
  for mrna in mRNA_species_array
    buffer *= "\n\tbackgroundGeneExpressionControlTermDict(\'$mrna\') = 0.1;  % $mrna"
  end
  # translation specific correction factor
  buffer *= "\n\n\t% translation specific correction factor"
  buffer *= "\n\ttranslationSpecificCorrectionFactor = containers.Map();"
  for pro in protein_species_array
    buffer *= "\n\ttranslationSpecificCorrectionFactor(\'$pro\') = 1.0;  % $pro"
  end
  # cooperativity number in transfer function
  buffer *= "\n\n\t% cooperativity number in transfer function"
  buffer *= "\n\tcooperativity = 1;"
  #buffer *= "\n\n\tbackground_mRNA_synthesis_rate_vector = 0.01*ones($length_TXTL)"
  # load txtl constants buffer
  buffer *= "\n"
  path_head = Base.@__DIR__
  if host_type == "bacteria"
    buffer *= replace(replace(m_include_constants_from_literature(joinpath(path_head,
        "txtl_constants_ecoli.jl"),"\n\t"), "#" => "%"), r"(  | \t|\t\t)% " => "  ; % ")
  else
    buffer *= replace(m_include_constants_from_literature(joinpath(path_head,
        "txtl_constants_hl60.jl"),"\n"), "#" => "%")
  end

  #---------------------------------
  # put all stuff in a Dictionary
  buffer *= "\n\n\t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
  buffer *= "\n\t% put all stuff in a Dictionary"
  buffer *= "\n\tdataDictionary = containers.Map();"
  buffer *= "\n\tdataDictionary(\'stoichiometric_matrix\') = stoichiometric_matrix;"
  buffer *= "\n\tdataDictionary(\'all_species_dict\') = all_species_dict;"
  buffer *= "\n\tdataDictionary(\'all_species_reversed_dict\') = all_species_reversed_dict;"
  buffer *= "\n\tdataDictionary(\'rnx_species_array\') = rnx_species_array;"
  buffer *= "\n\tdataDictionary(\'mRNA_species_array\') = mRNA_species_array;"
  buffer *= "\n\tdataDictionary(\'protein_species_array\') = protein_species_array;"
  buffer *= "\n\tdataDictionary(\'initial_condition\') = initial_condition;"

  buffer *= "\n\tdataDictionary(\'kcat_signaling\') = kcat_signaling;"
  buffer *= "\n\tdataDictionary(\'MonodAffinityConstantDict\') = MonodAffinityConstantDict;"
  buffer *= "\n\tdataDictionary(\'W_value_dict\') = W_value_dict;"
  buffer *= "\n\tdataDictionary(\'transferFunctionDisassociationConstantDict\') = transferFunctionDisassociationConstantDict;"
  buffer *= "\n\tdataDictionary(\'transcriptionSpecificCorrectionFactorDict\') = transcriptionSpecificCorrectionFactorDict;"
  buffer *= "\n\tdataDictionary(\'backgroundGeneExpressionControlTermDict\') = backgroundGeneExpressionControlTermDict;"
  buffer *= "\n\tdataDictionary(\'translationSpecificCorrectionFactor\') = translationSpecificCorrectionFactor;"
  buffer *= "\n\tdataDictionary(\'cooperativity\') = cooperativity;"
  buffer *= "\n\tdataDictionary(\'RNAPConcentration\') = rnapII_concentration;"
  buffer *= "\n\tdataDictionary(\'RIBOConcentration\') = ribosome_concentration;"
  buffer *= "\n\tdataDictionary(\'degradationConstantmRNA\') = degradation_constant_mRNA;"
  buffer *= "\n\tdataDictionary(\'degradationConstantProtein\') = degradation_constant_protein;"
  buffer *= "\n\tdataDictionary(\'kcatTranscription\') = kcat_transcription;"
  buffer *= "\n\tdataDictionary(\'kcatTranslation\') = kcat_translation;"
  buffer *= "\n\tdataDictionary(\'avgGeneConcentration\') = avg_gene_concentration;"
  buffer *= "\n\tdataDictionary(\'transcriptionSaturationConstant\') = saturation_transcription;"
  buffer *= "\n\tdataDictionary(\'translationSaturationConstant\') = saturation_translation;"
  buffer *= "\n\tdataDictionary(\'specificGrowthRate\') = maximum_specific_growth_rate - death_rate_constant;"

  # buffer *= "\n\n\treturn dataDictionary"
  buffer *= "\n\nend"

  return buffer
end

function m_build_simulation_buffer(NoExtracellularSpecies::Int64)
  # buffer = build_copyright_header_buffer()
  buffer = "\n% set up ODE, get the derivatives
function dydt = Balances(t,y,dataDictionary)
  % correct for negatives
  y(find(y < 0)) = 0.0;

  % load data from data dictionary
  rnx_species = dataDictionary(\'rnx_species_array\');
  mRNA_species = dataDictionary(\'mRNA_species_array\');
  protein_species = dataDictionary(\'protein_species_array\');
  stoichiometry = dataDictionary(\'stoichiometric_matrix\');
  all_species_dict = dataDictionary(\'all_species_dict\');  % String --> number
  specific_growth_rate = dataDictionary(\'specificGrowthRate\');
  transcription_correction = dataDictionary(\'transcriptionSpecificCorrectionFactorDict\');
  translation_correction = dataDictionary(\'translationSpecificCorrectionFactor\');
  degradation_constant_mRNA = dataDictionary(\'degradationConstantmRNA\');
  degradation_constant_protein = dataDictionary(\'degradationConstantProtein\');

  % kinetics
  [rnx_rate, transcription_rate, translation_rate] = Kinetics(y, dataDictionary);

  % calculate the derivatives
  dydt = zeros(1, length(all_species_dict));

  % signaling reaction network
  X2 = zeros(1, (length(rnx_species)-$(NoExtracellularSpecies)));  % species inside the cell
  for i = $(NoExtracellularSpecies+1):length(rnx_species)  % get X2 value from y
    X2(i - $(NoExtracellularSpecies)) = y(all_species_dict(char(rnx_species(i))));
  end
  dX1 = stoichiometry(1:$(NoExtracellularSpecies), :) * rnx_rate;  % extracellular metabolites
  dX2 = stoichiometry($(NoExtracellularSpecies)+1:end, :)*rnx_rate - transpose(specific_growth_rate*X2);  % INTRA
  dX = [dX1; dX2];
  for i = 1:length(rnx_species)  % return derivatives back to dydt
    dydt(all_species_dict(char(rnx_species(i)))) = dX(i);
  end
  dydt(1) = specific_growth_rate * y(1);  % for Biomass

  % TXTL network
  for i = 1:length(mRNA_species)
  % seems easy to add additional terms to account for the signaling effect
    mRNA_id = all_species_dict(char(mRNA_species(i)));
    p_id = all_species_dict(char(protein_species(i)));
    dydt(mRNA_id) = transcription_rate(char(mRNA_species(i))) - (specific_growth_rate + transcription_correction(char(mRNA_species(i)))*degradation_constant_mRNA) * y(mRNA_id);
    dydt(p_id) = translation_rate(char(protein_species(i))) - (specific_growth_rate + translation_correction(char(protein_species(i)))*degradation_constant_protein) * y(p_id);
  end
  dydt = dydt';
end"
  return buffer
end


function m_build_solveODEBalances_buffer(all_species_array::Array,
  all_species_dict::Dict{String, Int}, mRNA_species::Array, protein_species::Array)

  buffer = "function [t, y] = SolveBalances(TStart, TStop, TStep, dataDictionary)" *
    "\n\tt = TStart:TStep:TStop;"
  buffer *= "\n\tall_species_reversed_dict = dataDictionary(\'all_species_reversed_dict\');"
  # initial conditions
  buffer *= "\n\t% initial species concentration"
  buffer *= "\n\ty0 = dataDictionary(\'initial_condition\');"
  # run simulation
  buffer *= "\n\n\t[t, y] = ode23s(@(t,y) Balances(t,y,dataDictionary), t, y0);"
  # plotting
  buffer *= "\n\n\t% plotting results" *
    "\n\tsubplt_col = ceil(length(y0)/4);" *
    "\n\tfor i = 1:length(y0)" *
    "\n\t\tsubplot(4, subplt_col, i)" *
    "\n\t\tplot(t, y(:, i))" *
    "\n\t\ttitle(char(all_species_reversed_dict(i)), 'Interpreter', 'none')" *
    "\n\tend"
  buffer *= "\nend"

  return buffer
end


# FBA data dictionary generation
function m_generate_FBA_data_dictionary(all_rnx_list::Array,
  rnx_species_array::Array, extra_species_num::Int)
  secrete_id = Set{Int64}()  # For initialize the coefficient array
  # data dictionary
  buffer = "function data_dictionary = DataDictionary(time_start,time_stop,time_step)"
  # stoichiometry
  buffer *= "\n\n\t% load stoichiometry"
  buffer *= "\n\tstoichiometric_matrix = load(\'stoichiometry.dat\');"
  # Setup default flux bounds array
  buffer *= "\n\t% Setup default flux bounds array"
  buffer *= "\n\tdefault_bounds_array = ["
  for (id, rnx) in enumerate(all_rnx_list)
    buffer *= "\n\t\t0 100.0; % $(id) $(rnx.rnxName)"
  end
  buffer *= "\n\t];"
  # Setup default species bounds array
  buffer *= "\n\n\t% Setup default species bounds array"
  buffer *= "\n\tspecies_bounds_array = ["
  for (id, sp) in enumerate(rnx_species_array)
    if id > extra_species_num
      buffer *= "\n\t\t0.0 0.0;  % $(id) $(sp)"
    else
      buffer *= "\n\t\t-1.0 1.0; % $(id) $(sp)"
    end
  end
  buffer *= "\n\t];"
  # Setup the objective coefficient array
  buffer *= "\n\n\t% Setup the objective coefficient array"
  buffer *= "\n\tobjective_coefficient_array = ["
  for (id, rnx) in enumerate(all_rnx_list)
    buffer *= "\n\t\t0.0, ... % $(id) $(rnx.rnxName)"
    if rnx.rnxType == "secrete"
        push!(secrete_id, id)
    end
  end
  buffer *= "\n\t];"
  # List of reation strings - used to write flux report
  buffer *= "\n\n\t% List of reation strings - used to write flux report"
  buffer *= "\n\tlist_of_reaction_strings = {..."
  for rnx in all_rnx_list
    buffer *= "\n\t\t\'$(rnx.rnxName)\', ..."
  end
  buffer *= "\n\t};"
  # List of metabolite strings - used to write flux report
  buffer *= "\n\n\t% List of metabolite strings - used to write flux report"
  buffer *= "\n\tlist_of_metabolite_symbols = [..."
  for sp in rnx_species_array
    buffer *= "\n\t\t\'$(sp)\'..."
  end
  buffer *= "\n\t];"
  # return block -
  buffer *= "\n\n"
  buffer *= "\t% =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
  buffer *= "\tdata_dictionary = containers.Map();\n"
  buffer *= "\tdata_dictionary(\'stoichiometric_matrix\') = stoichiometric_matrix;\n"
  buffer *= "\tdata_dictionary(\'objective_coefficient_array\') = objective_coefficient_array;\n"
  buffer *= "\tdata_dictionary(\'default_flux_bounds_array\') = default_bounds_array;\n"
  buffer *= "\tdata_dictionary(\'species_bounds_array\') = species_bounds_array;\n"
  buffer *= "\tdata_dictionary(\'list_of_reaction_strings\') = list_of_reaction_strings;\n"
  buffer *= "\tdata_dictionary(\'list_of_metabolite_symbols\') = list_of_metabolite_symbols;\n"
  buffer *= "\tdata_dictionary(\'extra_species_num\') = $(extra_species_num);\n"
  buffer *= "\t% =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
  buffer *= "end\n"

  # maximize_product_dictionary
  buffer2 = "% maximize_product_dictionary"
  buffer2 *= "\nfunction dataDict = maximizeProductDictionary(time_start,time_stop,time_step)"
  # load the data dictionary -
  buffer2 *= "\n\t% load the data dictionary -"
  buffer2 *= "\n\tdataDict = DataDictionary(time_start,time_stop,time_step);"
  # Modify the data dictionary -
  buffer2 *= "\n\t% Modify the data dictionary -"
  buffer2 *= "\n\tobjective_coefficient_array = dataDict(\'objective_coefficient_array\');"
  secrete_id = sort(collect(secrete_id))
  for id in secrete_id
      buffer2 *= "\n\tobjective_coefficient_array($(id)) = -1;"
  end
  buffer2 *= "\n\tdataDict(\'objective_coefficient_array\') = objective_coefficient_array;"
  buffer2 *= "\nend"

  return buffer, buffer2
end


# FBA/FVA Solve.m
function m_build_FBA_solve_buffer()
  buffer =
"clc;
% load the data dictionary -
data_dictionary = maximizeProductDictionary(0,10,1);

% solve the lp problem -
[Flux,fVal,UptakeRate,EXITFLAG] = FluxDriver(data_dictionary, 1);
calculated_flux_array = FVA(data_dictionary, 1);

fprintf(\'objective value is %d \\n\', fVal);
fprintf(\'flux array: \\n\');
disp(Flux);
fprintf(\'uptake rate: \\n\');
disp(UptakeRate);
fprintf(\'FVA results: \\n\');
disp(calculated_flux_array);"
  return buffer
end


# rFBA rRules.m
function m_build_rRules_buffer()
  buffer =
"function new_bounds = rRules(bounds, varargin)
  % bounds: lower/upper bounds on flux
  % varargin: variable number of arguments, optional

  % initialize regulatory factors
  rfs = ones(size(bounds));
  % fill in regulatory rules below

  % update bounds
  new_bounds = rfs.*bounds;

end"
  return buffer
end


# rFBA Solve.m
function m_build_rFBA_solve_buffer()
  buffer =
"clc;
% load the data dictionary -
data_dictionary = maximizeProductDictionary(0,10,1);

% apply regulatory rules
new_bounds = rRules(data_dictionary(\'default_flux_bounds_array\'));
data_dictionary(\'default_flux_bounds_array\') = new_bounds;

% solve the lp problem -
[Flux,fVal,UptakeRate,EXITFLAG] = FluxDriver(data_dictionary, 1);
calculated_flux_array = FVA(data_dictionary, 1);

fprintf(\'objective value is %d \\n\', fVal);
fprintf(\'flux array: \\n\');
disp(Flux);
fprintf(\'uptake rate: \\n\');
disp(UptakeRate);
fprintf(\'FVA results: \\n\');
disp(calculated_flux_array);"
  return buffer
end
