# specific to Julia, that is, related to the output file format (language).

#=
F: read the file containing constants
I:
O:
BA
=#
function py3_include_constants_from_literature(src_file_name, pad_string)
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
function py3_build_kinetics_buffer(all_species_dict::Dict{String, Int},
  all_rnx_list::Array, all_txtl_dict::Dict, sys2user::Dict)
  kinetics = "import numpy as np" *
             "\n\ndef calculate_kinetics(X, data_dictionary):" *
             "\n\n\t# load all species dictionary" *
             "\n\tall_species_dict = data_dictionary[\"all_species_dict\"]  # String-->Int"
  # for signaling
  kinetics *= "\n\n\t# generate kinetics for signaling, assume saturation kinetics"
  totalRnxNo = length(all_rnx_list)
  kinetics *= "\n\t# load Monod affinity constants" *
              "\n\tMonodK = data_dictionary[\"MonodAffinityConstantDict\"]  # String-->Float"
  MonodAffinityConstant_String_Array = Array{String,1}()  # collect monod affinity constant symbol
  kinetics *= "\n\t# load reaction kinetic constants" *
              "\n\tkcat = data_dictionary[\"kcat_signaling\"]  # in order of rnx" *
              "\n\n\t# write reaction rate equations" *
              "\n\trnx_rate_vector = np.zeros($totalRnxNo)"
  for (index, rnx) in enumerate(all_rnx_list)  # go thru every rnx
    tmp_rnx_rate = "kcat[$(index - 1)]"
    if isdefined(rnx, :catalysts)  # enzyme
      for token in rnx.catalysts
        tmp_rnx_rate *= "*X[$(all_species_dict[token.oriBioName] - 1)]"
      end
    end
    if isdefined(rnx,:reactants)  # substrate
      for token in rnx.reactants
        tmp_rnx_rate *= "*X[$(all_species_dict[token.oriBioName] - 1)]/(MonodK[\"MonodK~rnx$(index)~$(token.oriBioName)\"]+X[$(all_species_dict[token.oriBioName] - 1)])"
        push!(MonodAffinityConstant_String_Array, "MonodK~rnx$(index)~$(token.oriBioName)")
      end
    end
    kinetics *= "\n\trnx_rate_vector[$(index - 1)] = "*tmp_rnx_rate
  end


  # for TXTL
  W_string_array = Array{String,1}()  # W: weight of protein action on transcription
  disassociation_const_string_array = Array{String,1}()  # disassociation constant in tranfer function
  kinetics *= "\n\n\t# generate kinetics & control terms for TXTL" *
              "\n\t# load data from data dictionary" *
              "\n\tW_value_dict = data_dictionary[\"W_value_dict\"]  #String-->Float" *
              "\n\tcoop = data_dictionary[\"cooperativity\"]" *
              "\n\tdisassociation_constant_dict = data_dictionary[\"transferFunctionDisassociationConstantDict\"]" *
              "\n\tbackground_control_term = data_dictionary[\"backgroundGeneExpressionControlTermDict\"]" *
              "\n\tkcatTX = data_dictionary[\"kcatTranscription\"]" *
              "\n\tRNAPconc = data_dictionary[\"RNAPConcentration\"]" *
              "\n\tgeneConc = data_dictionary[\"avgGeneConcentration\"]" *
              "\n\tsaturationTX = data_dictionary[\"transcriptionSaturationConstant\"]" *
              "\n\tkcatTL = data_dictionary[\"kcatTranslation\"]" *
              "\n\tRIBOconc = data_dictionary[\"RIBOConcentration\"]" *
              "\n\tsaturationTL = data_dictionary[\"translationSaturationConstant\"]"
  totalTXTLNo = length(all_txtl_dict)
  kinetics *= "\n\n\t# write control terms, transcription rate equations, and translation rate equations" *
             "\n\tTX_control_term = {}  # control term: string -> float" *
             "\n\tTX_rate_vector = {}  # transcription rate: string -> float" *
             "\n\tTL_rate_vector = {}  # translation rate: string -> float"
  for (key, txtl) in all_txtl_dict # go thru every txtl
    tmp_protein_string = replace(key, sys2user["MRNA"] => sys2user["PROTEIN"], count=1)
    kinetics *= "\n\t# $key and $tmp_protein_string"
    up_factors_array = Array{String,1}()  # collection of upregulation factors name: targetedmRNA_factor(s)
    if !isempty(txtl.activationProtein)  # upregulation
      for token_array in txtl.activationProtein
        tmp_up = ""
        tmp_up_name = key*"_"
        for token in token_array
          tmp_W = "W~$(token.oriBioName)~$(key)"
          tmp_up *= "W_value_dict[\"$(tmp_W)\"] * X[$(all_species_dict[token.oriBioName] - 1)]**coop / ((disassociation_constant_dict[\"KD~$(token.oriBioName)\"])**coop+X[$(all_species_dict[token.oriBioName] - 1)]**coop)*"
          push!(W_string_array, tmp_W)
          push!(disassociation_const_string_array, "KD~$(token.oriBioName)")
          tmp_up_name *= token.oriBioName*"_"
        end
        tmp_up = chop(tmp_up)
        tmp_up_name = chop(tmp_up_name)
        push!(up_factors_array, tmp_up_name)
        kinetics *= "\n\t$(tmp_up_name) = $(tmp_up)"
      end
    end
    down_factors_array = Array{String,1}()  # collection of downregulation factors name: targetedmRNA_factor(s)
    if !isempty(txtl.inhibitionProtein)  # downregulation
      for token_array in txtl.inhibitionProtein
        tmp_down = ""
        tmp_down_name = key*"_"
        for token in token_array
          tmp_W = "W~$(token.oriBioName)~$(key)"
          tmp_down *= "W_value_dict[\"$(tmp_W)\"] * X[$(all_species_dict[token.oriBioName] - 1)]**coop/((disassociation_constant_dict[\"KD~$(token.oriBioName)\"])**coop+X[$(all_species_dict[token.oriBioName] - 1)]**coop)*"
          push!(W_string_array, tmp_W)
          push!(disassociation_const_string_array, "KD~$(token.oriBioName)")
          tmp_down_name *= token.oriBioName*"_"
        end
        tmp_down = chop(tmp_down)
        tmp_down_name = chop(tmp_down_name)
        push!(down_factors_array, tmp_down_name)
        kinetics *= "\n\t$(tmp_down_name) = $(tmp_down)"
      end
    end
    # combine to form up and down, then merge up and down to form control term
    if !isempty(up_factors_array) && !isempty(down_factors_array)
      up_action = ""
      for up in up_factors_array
        up_action *= up*" +"
      end
      kinetics *= "\n\tup_action_on_$key = $(chop(up_action))"
      down_action = ""
      for down in down_factors_array
        down_action *= down*" +"
      end
      kinetics *= "\n\tdown_action_on_$key = $(chop(down_action))" *
                  "\n\tTX_control_term[\"$key\"] = (background_control_term[\"$key\"] + up_action_on_$key)/(1 + background_control_term[\"$key\"] + up_action_on_$key + down_action_on_$key)"
    elseif !isempty(up_factors_array)
      up_action = ""
      for up in up_factors_array
        up_action *= up*" +"
      end
      kinetics *= "\n\tup_action_on_$key = $(chop(up_action))" *
                  "\n\tTX_control_term[\"$key\"] = (background_control_term[\"$key\"] + up_action_on_$key)/(1 + background_control_term[\"$key\"] + up_action_on_$key)"
    elseif !isempty(down_factors_array)
      down_action = ""
      for down in down_factors_array
        down_action *= down*" +"
      end
      kinetics *= "\n\tdown_action_on_$key = $(chop(down_action))" *
                  "\n\tTX_control_term[\"$key\"] = background_control_term[\"$key\"]/(1 + background_control_term[\"$key\"] + down_action_on_$key)"
    else
      # seems impossible, a TXTL dict is constructed b/c something acts on the key
    end
    kinetics *= "\n\tTX_rate_vector[\"$key\"] = TX_control_term[\"$key\"] * kcatTX * RNAPconc*geneConc/(saturationTX+geneConc)" *
                "\n\tTL_rate_vector[\"$tmp_protein_string\"] = kcatTL * RIBOconc*X[$(all_species_dict[key] - 1)]/(saturationTL+X[$(all_species_dict[key] - 1)])"
  end

  kinetics *= "\n\n\treturn rnx_rate_vector, TX_rate_vector, TL_rate_vector"

  return (kinetics, MonodAffinityConstant_String_Array, W_string_array, disassociation_const_string_array)
end

function py3_build_data_dictionary_buffer(host_type::AbstractString, all_species_array::Array,
  all_species2index_dict::Dict,
  rnx_species_array::Array, all_rnx_list::Array, all_txtl_dict::Dict,
  Monod_affinity_constant_array::Array, W_string_array::Array,
  disassociation_const_array::Array, mRNA_species_array::Array, protein_species_array::Array)
  # data dictionary
  buffer = "import numpy as np" *
    "\nimport math"
  buffer *= "\n\ndef generate_model_parameters_dictionary():"
  # stoichiometry
  buffer *= "\n\n\t# load stoichiometry" *
    "\n\twith open(\"./stoichiometry.dat\", \"r\") as file:" *
    "\n\t\tmatrix = [[float(s) for s in line.strip().split()] for line in file]" *
    "\n\t\tprint(np.shape(matrix))" *
    "\n\tstoichiometric_matrix = np.matrix(matrix)" *
    "\n\tprint(np.shape(matrix))"

  # initial condition
  No_species = length(all_species_array)
  # initial concentration
  buffer *= "\n\n\t# initial species concentration" *
    "\n\tinitial_condition = [\n\t\t"
  buffer *= join(["0.001,  # $(i-1)  $(all_species_array[i])" for i = 1:(No_species - 1)], "\n\t\t")
  buffer *= "\n\t\t0.001   # $(No_species-1)  $(all_species_array[No_species])"
  buffer *= "\n\t]"

  # all species dictionary
  buffer *= "\n\n\t# all species dictionary"
  buffer *= "\n\tall_species_dict = {}  # String --> Int"
  #all_species_dict = Dict{String, Int64}()  # create all species dict: string --> int64
  for (st, index) in all_species2index_dict
    buffer *= "\n\tall_species_dict[\"$st\"] = $(index - 1)"
    #all_species_dict[st] = index
  end
  # all species dictionary_reversed
  buffer *= "\n\n\t# all species dictionary_reversed"
  buffer *= "\n\tall_species_reversed_dict = {}  # Int --> String"
  buffer *= "\n\tfor key, val in all_species_dict.items():"
  buffer *= "\n\t\tall_species_reversed_dict[val] = key"
  # all rnx species array
  buffer *= "\n\n\t# all rnx species array"
  buffer *= "\n\trnx_species_array = ["
  for st in rnx_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"
  # all mRNA species array
  buffer *= "\n\n\t# all mRNA species array"
  buffer *= "\n\tmRNA_species_array = ["
  for st in mRNA_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"
  # all protein species array
  buffer *= "\n\n\t# all protein species array"
  buffer *= "\n\tprotein_species_array = ["
  for st in protein_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"

  # signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)
  buffer *= "\n\n\t# signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)"
  buffer *= "\n\tkcat_signaling = 1e-3*np.ones($(length(all_rnx_list)))  # kcat[#reaction]: reaction name"
  for (id, rnx) in enumerate(all_rnx_list)
    buffer*= "\n\tkcat_signaling[$(id - 1)] = 1.0  # kcat: $(rnx.rnxName)"
  end
  # Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate
  buffer *= "\n\n\t# Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate"
  buffer *= "\n\tMonodAffinityConstantDict = {}"
  for MK in Monod_affinity_constant_array
    buffer *= "\n\tMonodAffinityConstantDict[\"$MK\"] = 1.0"
  end

  # W_value in transcription control term, nomenclature: W_targetmRNA_actor
  buffer *= "\n\n\n\t# W_value in transcription control term, nomenclature: W_targetmRNA_actor"
  buffer *= "\n\tW_value_dict = {}"
  for WS in W_string_array
    buffer *= "\n\tW_value_dict[\"$WS\"] = 1.0"
  end
  # disassociation constants in transfer function, nomenclature: KD_speciesName
  buffer *= "\n\n\t# disassociation constants in transfer function, nomenclature: KD_speciesName"
  buffer *= "\n\ttransferFunctionDisassociationConstantDict = {}"

  for DC in disassociation_const_array
    buffer *= "\n\ttransferFunctionDisassociationConstantDict[\"$DC\"] = 1.0"
  end
  # transcription specific correction factor
  buffer *= "\n\n\t# transcription specific correction factor"
  buffer *= "\n\ttranscriptionSpecificCorrectionFactorDict = {}"
  for mrna in mRNA_species_array
    buffer *= "\n\ttranscriptionSpecificCorrectionFactorDict[\"$mrna\"] = 1.0  # $mrna"
  end
  # background gene expression control term
  buffer *= "\n\n\t# background gene expression control term"
  buffer *= "\n\tbackgroundGeneExpressionControlTermDict = {}"
  for mrna in mRNA_species_array
    buffer *= "\n\tbackgroundGeneExpressionControlTermDict[\"$mrna\"] = 0.1  # $mrna"
  end
  # translation specific correction factor
  buffer *= "\n\n\t# translation specific correction factor"
  buffer *= "\n\ttranslationSpecificCorrectionFactor = {}"
  for pro in protein_species_array
    buffer *= "\n\ttranslationSpecificCorrectionFactor[\"$pro\"] = 1.0  # $pro"
  end
  # cooperativity number in transfer function
  buffer *= "\n\n\t# cooperativity number in transfer function"
  buffer *= "\n\tcooperativity = 1"
  #buffer *= "\n\n\tbackground_mRNA_synthesis_rate_vector = 0.01*ones($length_TXTL)"
  # load txtl constants buffer
  buffer *= "\n\n \n"
  if host_type == "bacteria"
    buffer *= replace(replace(py3_include_constants_from_literature(
      joinpath(Base.@__DIR__, "txtl_constants_ecoli.jl"),"\n\t"),
      ")^3" => ")**3"), "log(" => "math.log(")
  else
    buffer *= replace(replace(py3_include_constants_from_literature(
      joinpath(Base.@__DIR__, "txtl_constants_hl60.jl"),"\n\t"),
      ")^3" => ")**3"), "log(" => "math.log(")
  end

  #---------------------------------
  # put all stuff in a Dictionary
  buffer *= "\n\n\t####################################"
  buffer *= "\n\t# put all stuff in a Dictionary"
  buffer *= "\n\tdataDictionary = {}  # string -> any"
  buffer *= "\n\tdataDictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix"
  buffer *= "\n\tdataDictionary[\"all_species_dict\"] = all_species_dict"
  buffer *= "\n\tdataDictionary[\"all_species_reversed_dict\"] = all_species_reversed_dict"
  buffer *= "\n\tdataDictionary[\"rnx_species_array\"] = rnx_species_array"
  buffer *= "\n\tdataDictionary[\"mRNA_species_array\"] = mRNA_species_array"
  buffer *= "\n\tdataDictionary[\"protein_species_array\"] = protein_species_array"
  buffer *= "\n\tdataDictionary[\"initial_condition\"] = initial_condition"

  buffer *= "\n\tdataDictionary[\"kcat_signaling\"] = kcat_signaling"
  buffer *= "\n\tdataDictionary[\"MonodAffinityConstantDict\"] = MonodAffinityConstantDict"
  buffer *= "\n\tdataDictionary[\"W_value_dict\"] = W_value_dict"
  buffer *= "\n\tdataDictionary[\"transferFunctionDisassociationConstantDict\"] = transferFunctionDisassociationConstantDict"
  buffer *= "\n\tdataDictionary[\"transcriptionSpecificCorrectionFactorDict\"] = transcriptionSpecificCorrectionFactorDict"
  buffer *= "\n\tdataDictionary[\"backgroundGeneExpressionControlTermDict\"] = backgroundGeneExpressionControlTermDict"
  buffer *= "\n\tdataDictionary[\"translationSpecificCorrectionFactor\"] = translationSpecificCorrectionFactor"
  buffer *= "\n\tdataDictionary[\"cooperativity\"] = cooperativity"
  buffer *= "\n\tdataDictionary[\"RNAPConcentration\"] = rnapII_concentration"
  buffer *= "\n\tdataDictionary[\"RIBOConcentration\"] = ribosome_concentration"
  buffer *= "\n\tdataDictionary[\"degradationConstantmRNA\"] = degradation_constant_mRNA"
  buffer *= "\n\tdataDictionary[\"degradationConstantProtein\"] = degradation_constant_protein"
  buffer *= "\n\tdataDictionary[\"kcatTranscription\"] = kcat_transcription"
  buffer *= "\n\tdataDictionary[\"kcatTranslation\"] = kcat_translation"
  buffer *= "\n\tdataDictionary[\"avgGeneConcentration\"] = avg_gene_concentration"
  buffer *= "\n\tdataDictionary[\"transcriptionSaturationConstant\"] = saturation_transcription"
  buffer *= "\n\tdataDictionary[\"translationSaturationConstant\"] = saturation_translation"
  buffer *= "\n\tdataDictionary[\"specificGrowthRate\"] = maximum_specific_growth_rate - death_rate_constant"

  buffer *= "\n\n\treturn dataDictionary"

  return buffer
end

function py3_build_simulation_buffer(NoExtracellularSpecies::Int64)
  # buffer = build_copyright_header_buffer()
  buffer = "from Kinetics import calculate_kinetics" *
    "\nimport numpy as np" *
  "\n# set up ODE, get the derivatives
def Balances(y, t, dataDictionary):
  # correct for negatives
  y = [x if x > 0 else 0 for x in y]

  # load data from data dictionary
  rnx_species = dataDictionary[\"rnx_species_array\"]
  mRNA_species = dataDictionary[\"mRNA_species_array\"]
  protein_species = dataDictionary[\"protein_species_array\"]
  stoichiometry = dataDictionary[\"stoichiometric_matrix\"]
  all_species_dict = dataDictionary[\"all_species_dict\"]  # String --> number
  specific_growth_rate = dataDictionary[\"specificGrowthRate\"]
  transcription_correction = dataDictionary[\"transcriptionSpecificCorrectionFactorDict\"]
  translation_correction = dataDictionary[\"translationSpecificCorrectionFactor\"]
  degradation_constant_mRNA = dataDictionary[\"degradationConstantmRNA\"]
  degradation_constant_protein = dataDictionary[\"degradationConstantProtein\"]

  # kinetics
  rnx_rate, transcription_rate, translation_rate = calculate_kinetics(y, dataDictionary)

  # calculate the derivatives
  dydt = np.zeros(len(all_species_dict))

  # signaling reaction network
  X2 = np.zeros(len(rnx_species)-$(NoExtracellularSpecies))  # species inside the cell
  for id, st in enumerate(rnx_species[$(NoExtracellularSpecies):]):
    # get X2 value from y
    X2[id] = y[all_species_dict[st]]

  dX1 = np.dot(stoichiometry[0:$(NoExtracellularSpecies), :], rnx_rate)  # extracellular metabolites
  dX2 = np.dot(stoichiometry[$(NoExtracellularSpecies): , :], rnx_rate) - specific_growth_rate * X2  # INTRA
  # print(np.shape(dX1), np.shape(dX2))
  dX = np.concatenate((dX1, dX2), axis=1)
  # print(np.shape(dX))
  for id, st in enumerate(rnx_species):
    # return derivatives back to dydt
    dydt[all_species_dict[st]] = dX[0, id]

  dydt[0] = specific_growth_rate * y[0]  # for Biomass

  # TXTL network
  for i in range(len(mRNA_species)):
  # seems easy to add additional terms to account for the signaling effect
    mRNA_id = all_species_dict[mRNA_species[i]]
    p_id = all_species_dict[protein_species[i]]
    dydt[mRNA_id] = transcription_rate[mRNA_species[i]] - (specific_growth_rate + transcription_correction[mRNA_species[i]]*degradation_constant_mRNA) * y[mRNA_id]
    dydt[p_id] = translation_rate[protein_species[i]] - (specific_growth_rate + translation_correction[protein_species[i]]*degradation_constant_protein) * y[p_id]


  return dydt"

  return buffer
end


function py3_build_solveODEBalances_buffer(all_species_array::Array,
  all_species_dict::Dict{String, Int}, mRNA_species::Array, protein_species::Array)

  buffer = "import scipy.integrate as integrate" *
    "\nimport matplotlib.pyplot as plt" *
    "\nimport numpy as np" *
    "\nimport math" *
    "\nfrom DataDictionary import generate_model_parameters_dictionary" *
    "\nfrom Kinetics import calculate_kinetics" *
    "\nfrom Balances import Balances" *
    "\n" *
    "\ndef SolveBalances(TStart, TStop, TStep, dataDictionary):" *
    "\n\ttimeNum = int((TStop - TStart)/TStep)" *
    "\n\tt = np.linspace(TStart, TStop, timeNum)"

  # data dataDictionary
  buffer *= "\n\tall_species_reversed_dict = dataDictionary[\"all_species_reversed_dict\"]"
  buffer *= "\n\tinitial_condition = dataDictionary[\"initial_condition\"]"
  buffer *= "\n\ty0 = dataDictionary[\"initial_condition\"]"
  # run simulation
  buffer *= "\n\t# Calling the solver"
  buffer *= "\n\ty = integrate.odeint(Balances, y0, t, args=(dataDictionary,))"
  # plotting
  buffer *= "\n\n\t# plotting results" *
    "\n\tcol = len(y0)" *
    "\n\tsubplt_col = math.ceil(col/4)" *
    "\n\tfor i in range(col):" *
    "\n\t\tplt.subplot(4, subplt_col, i+1)" *
    "\n\t\tplt.plot(t, y[:, i])" *
    "\n\t\tplt.title(all_species_reversed_dict[i])" *
    "\n\tplt.show()"
  buffer *= "\n\treturn t, y"

  # add a test
  buffer *= "\n\n\n"
  buffer *= "## uncomment to test"
  buffer *= "\n# dataDictionary = generate_model_parameters_dictionary()"
  buffer *= "\n# SolveBalances(0, 10, 0.1, dataDictionary)"

  return buffer
end


# FBA data dictionary generation
function py3_generate_FBA_data_dictionary(all_rnx_list::Array,
  rnx_species_array::Array, extra_species_num::Int)
  secrete_id = Set{Int64}()  # For initialize the coefficient array
  # data dictionary
  buffer = "import numpy as np" *
    "\nimport math"
  buffer *= "\n\ndef DataDictionary(time_start,time_stop,time_step):"
  # stoichiometry
  buffer *= "\n\n\t# load stoichiometry" *
    "\n\twith open(\"./stoichiometry.dat\", \"r\") as file:" *
    "\n\t\tmatrix = [[float(s) for s in line.strip().split()] for line in file]" *
    # "\n\t\tprint(np.shape(matrix))" *
    "\n\tstoichiometric_matrix = np.matrix(matrix)"
    # "\n\tprint(np.shape(matrix))"

  # Setup default flux bounds array
  buffer *= "\n\t# Setup default flux bounds array"
  buffer *= "\n\tdefault_flux_bounds_array = ["
  for (id, rnx) in enumerate(all_rnx_list)
    buffer *= "\n\t\t[0, 100.0], # $(id) $(rnx.rnxName)"
  end
  buffer *= "\n\t]"
  # Setup default species bounds array
  buffer *= "\n\n\t# Setup default species bounds array"
  buffer *= "\n\tspecies_bounds_array = ["
  for (id, sp) in enumerate(rnx_species_array)
    if id > extra_species_num
      buffer *= "\n\t\t[0.0, 0.0], # $(id) $(sp)"
    else
      buffer *= "\n\t\t[-1.0, 1.0], # $(id) $(sp)"
    end
  end
  buffer *= "\n\t]"
  # Setup the objective coefficient array
  buffer *= "\n\n\t# Setup the objective coefficient array"
  buffer *= "\n\tobjective_coefficient_array = ["
  for (id, rnx) in enumerate(all_rnx_list)
    buffer *= "\n\t\t0.0, # $(id) $(rnx.rnxName)"
    if rnx.rnxType == "secrete"
        push!(secrete_id, (id-1))
    end
  end
  buffer *= "\n\t]"
  # List of reation strings - used to write flux report
  buffer *= "\n\n\t# List of reation strings - used to write flux report"
  buffer *= "\n\tlist_of_reaction_strings = ["
  for rnx in all_rnx_list
    buffer *= "\n\t\t\"$(rnx.rnxName)\", "
  end
  buffer *= "\n\t]"
  # List of metabolite strings - used to write flux report
  buffer *= "\n\n\t# List of metabolite strings - used to write flux report"
  buffer *= "\n\tlist_of_metabolite_symbols = ["
  for sp in rnx_species_array
    buffer *= "\n\t\t\"$(sp)\","
  end
  buffer *= "\n\t]"
  # return block -
  buffer *= "\n\n"
  buffer *= "\t# =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
  buffer *= "\tdata_dictionary = {}\n"
  buffer *= "\tdata_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix\n"
  buffer *= "\tdata_dictionary[\"objective_coefficient_array\"] = objective_coefficient_array\n"
  buffer *= "\tdata_dictionary[\"default_flux_bounds_array\"] = default_flux_bounds_array;\n"
  buffer *= "\tdata_dictionary[\"species_bounds_array\"] = species_bounds_array\n"
  buffer *= "\tdata_dictionary[\"list_of_reaction_strings\"] = list_of_reaction_strings\n"
  buffer *= "\tdata_dictionary[\"list_of_metabolite_symbols\"] = list_of_metabolite_symbols\n"
  buffer *= "\tdata_dictionary[\"extra_species_num\"] = $(extra_species_num)\n"
  buffer *= "\t# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
  buffer *= "\treturn data_dictionary\n"


  # maximize_product_dictionary
  buffer *= "\n\n\n# maximize_product_dictionary"
  buffer *= "\ndef maximize_product_dictionary(time_start,time_stop,time_step):"
  # load the data dictionary -
  buffer *= "\n\t# load the data dictionary -"
  buffer *= "\n\tdata_dictionary = DataDictionary(time_start,time_stop,time_step)"
  # Modify the data dictionary -
  buffer *= "\n\t# Modify the data dictionary -"
  buffer *= "\n\tobjective_coefficient_array = data_dictionary[\"objective_coefficient_array\"]"
  secrete_id = sort(collect(secrete_id))
  for id in secrete_id
      buffer *= "\n\tobjective_coefficient_array[$(id)] = -1"
  end
  buffer *= "\n\t# return -"
  buffer *= "\n\treturn data_dictionary"

  return buffer
end


# FBA/FVA Solve.py
function py3_build_FBA_solve_buffer()
  buffer =
"from FluxDriver import FluxDriver
from DataDictionary import maximize_product_dictionary
from FVA import calculate_flux_variabilty

# load the data dictionary -
data_dictionary = maximize_product_dictionary(0,10,1)

# solve the lp problem -
objective_value, flux_array, uptake_array, exit_flag = FluxDriver(data_dictionary)
calculated_flux_array = calculate_flux_variabilty(data_dictionary)

print(\"objective_value: \", objective_value)
print(\"flux array: \", flux_array)
print(\"uptake array: \", uptake_array)
print(\"FVA results: \\n\", calculated_flux_array)"
  return buffer
end


# rFBA rRules.jl
function py3_build_rRules_buffer()
  buffer = "import numpy as np
def rRules(bounds, *args):
  # bounds: lower/upper bounds on flux
  # args: variable number of arguments, optional
  bounds = np.array(bounds)

  # initialize regulatory factors
  rfs = np.ones(bounds.shape)
  # fill in regulatory rules below

  # update bounds
  new_bounds = rfs*bounds

  return new_bounds"
  return buffer
end


# rFBA Solve.py
function py3_build_rFBA_solve_buffer()
  buffer =
"from FluxDriver import FluxDriver
from DataDictionary import maximize_product_dictionary
from FVA import calculate_flux_variabilty
from rRules import rRules

# load the data dictionary -
data_dictionary = maximize_product_dictionary(0,10,1)

# apply regulatory rules
new_bounds = rRules(data_dictionary[\"default_flux_bounds_array\"])
data_dictionary[\"default_flux_bounds_array\"] = new_bounds

# solve the lp problem -
objective_value, flux_array, uptake_array, exit_flag = FluxDriver(data_dictionary)
calculated_flux_array = calculate_flux_variabilty(data_dictionary)

print(\"objective_value: \", objective_value)
print(\"flux array: \", flux_array)
print(\"uptake array: \", uptake_array)
print(\"FVA results: \\n\", calculated_flux_array)"
  return buffer
end
