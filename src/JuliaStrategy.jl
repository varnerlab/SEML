# specific to Julia, that is, related to the output file format (language).

#=
F: read the file containing constants
I:
O:
BA
=#
function jl_include_constants_from_literature(src_file_name, pad_string)
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
function jl_build_kinetics_buffer(all_species_dict::Dict{String, Int},
  all_rnx_list::Array, all_txtl_dict::Dict, sys2user::Dict)
  kinetics = "function calculate_kinetics(X, data_dictionary)" *
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
              "\n\trnx_rate_vector = zeros($totalRnxNo)"
  for (index, rnx) in enumerate(all_rnx_list)  # go thru every rnx
    tmp_rnx_rate = "kcat[$index]"
    if isdefined(rnx, :catalysts)  # enzyme
      for token in rnx.catalysts
        tmp_rnx_rate *= "*X[$(all_species_dict[token.oriBioName])]"
      end
    end
    if isdefined(rnx,:reactants)  # substrate
      for token in rnx.reactants
        tmp_rnx_rate *= "*X[$(all_species_dict[token.oriBioName])]/(MonodK[\"MonodK~rnx$(index)~$(token.oriBioName)\"]+X[$(all_species_dict[token.oriBioName])])"
        push!(MonodAffinityConstant_String_Array, "MonodK~rnx$(index)~$(token.oriBioName)")
      end
    end
    kinetics *= "\n\trnx_rate_vector[$index] = "*tmp_rnx_rate
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
             "\n\tTX_control_term = Dict{String, Float64}()  # control term" *
             "\n\tTX_rate_vector = Dict{String, Float64}()  # transcription rate" *
             "\n\tTL_rate_vector = Dict{String, Float64}()  # translation rate"
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
          tmp_up *= "W_value_dict[\"$(tmp_W)\"] * X[$(all_species_dict[token.oriBioName])]^coop/((disassociation_constant_dict[\"KD~$(token.oriBioName)\"])^coop+X[$(all_species_dict[token.oriBioName])]^coop)*"
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
          tmp_down *= "W_value_dict[\"$(tmp_W)\"] * X[$(all_species_dict[token.oriBioName])]^coop/((disassociation_constant_dict[\"KD~$(token.oriBioName)\"])^coop+X[$(all_species_dict[token.oriBioName])]^coop)*"
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
    kinetics *= "\n\tTX_rate_vector[\"$key\"] = TX_control_term[\"$key\"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)" *
                "\n\tTL_rate_vector[\"$tmp_protein_string\"] = kcatTL*RIBOconc*X[$(all_species_dict[key])]/(saturationTL+X[$(all_species_dict[key])])"
  end

  kinetics *= "\n\n\treturn (rnx_rate_vector, TX_rate_vector, TL_rate_vector)" *
              "\nend"
  return (kinetics, MonodAffinityConstant_String_Array, W_string_array, disassociation_const_string_array)
end

function jl_build_data_dictionary_buffer(host_type::AbstractString, all_species_array::Array,
  all_species2index_dict::Dict,
  rnx_species_array::Array, all_rnx_list::Array, all_txtl_dict::Dict,
  Monod_affinity_constant_array::Array, W_string_array::Array,
  disassociation_const_array::Array, mRNA_species_array::Array, protein_species_array::Array)
  # data dictionary
  buffer = "function generate_model_parameters_dictionary()"
  # stoichiometry
  buffer *= "\n\t# load stoichiometry"
  buffer *= "\n\tstoichiometric_matrix = readdlm(\"./stoichiometry.dat\")"

  No_species = length(all_species_array)
  # initial concentration
  buffer *= "\n\n\t# initial species concentration" *
    "\n\tinitial_condition = [\n\t\t"
  buffer *= join(["0.0,  # $(i)  $(all_species_array[i])" for i = 1:(No_species - 1)], "\n\t\t")
  buffer *= "\n\t\t0.0   # $(No_species) $(all_species_array[No_species])"
  buffer *= "\n\t]"

  # all species dictionary
  buffer *= "\n\n\t# all species dictionary"
  buffer *= "\n\tall_species_dict = Dict{String, Int64}()  # String --> Int"
  #all_species_dict = Dict{String, Int64}()  # create all species dict: string --> int64
  for (st, index) in all_species2index_dict
    buffer *= "\n\tall_species_dict[\"$st\"] = $index"
    #all_species_dict[st] = index
  end
  # all species dictionary_reversed
  buffer *= "\n\n\t# all species dictionary_reversed"
  buffer *= "\n\tall_species_reversed_dict = Dict{Int64, String}()  # Int --> String"
  buffer *= "\n\tfor (key,val) in all_species_dict"
  buffer *= "\n\t\tall_species_reversed_dict[val] = key"
  buffer *= "\n\tend"
  # all rnx species array
  buffer *= "\n\n\t# all rnx species array"
  buffer *= "\n\trnx_species_array::Array{String,1} = ["
  for st in rnx_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"
  # all mRNA species array
  buffer *= "\n\n\t# all mRNA species array"
  buffer *= "\n\tmRNA_species_array::Array{String,1} = ["
  for st in mRNA_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"
  # all protein species array
  buffer *= "\n\n\t# all protein species array"
  buffer *= "\n\tprotein_species_array::Array{String,1} = ["
  for st in protein_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"

  # signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)
  buffer *= "\n\n\t# signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)"
  buffer *= "\n\tkcat_signaling = ones($(length(all_rnx_list)))  # kcat[#reaction]: reaction name"
  for (id, rnx) in enumerate(all_rnx_list)
    buffer*= "\n\tkcat_signaling[$id] = 1e-4  # kcat: $(rnx.rnxName)"
  end
  # Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate
  buffer *= "\n\n\t# Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate"
  buffer *= "\n\tMonodAffinityConstantDict = Dict{String, Float64}()"
  for MK in Monod_affinity_constant_array
    # buffer *= "\n\t$MK = 1.0"
    # buffer *= "\n\tMonodAffinityConstant_dict[\"$MK\"] = $MK"
    buffer *= "\n\tMonodAffinityConstantDict[\"$MK\"] = 0.1"
  end

  # W_value in transcription control term, nomenclature: W_targetmRNA_actor
  buffer *= "\n\n\n\t# W_value in transcription control term, nomenclature: W_targetmRNA_actor"
  buffer *= "\n\tW_value_dict = Dict{String, Float64}()"
  for WS in W_string_array
    buffer *= "\n\tW_value_dict[\"$WS\"] = 1.0"
  end
  # disassociation constants in transfer function, nomenclature: KD_speciesName
  buffer *= "\n\n\t# disassociation constants in transfer function, nomenclature: KD_speciesName"
  buffer *= "\n\ttransferFunctionDisassociationConstantDict = Dict{String, Float64}()"
  for DC in disassociation_const_array
    buffer *= "\n\ttransferFunctionDisassociationConstantDict[\"$DC\"] = 1.0"
  end
  # transcription specific correction factor
  buffer *= "\n\n\t# transcription specific correction factor"
  buffer *= "\n\ttranscriptionSpecificCorrectionFactorDict = Dict{String, Float64}()"
  for mrna in mRNA_species_array
    buffer *= "\n\ttranscriptionSpecificCorrectionFactorDict[\"$mrna\"] = 1.0  # $mrna"
  end
  # background gene expression control term
  buffer *= "\n\n\t# background gene expression control term"
  buffer *= "\n\tbackgroundGeneExpressionControlTermDict = Dict{String, Float64}()"
  for mrna in mRNA_species_array
    buffer *= "\n\tbackgroundGeneExpressionControlTermDict[\"$mrna\"] = 0.001  # $mrna"
  end
  # translation specific correction factor
  buffer *= "\n\n\t# translation specific correction factor"
  buffer *= "\n\ttranslationSpecificCorrectionFactor = Dict{String, Float64}()"
  for pro in protein_species_array
    buffer *= "\n\ttranslationSpecificCorrectionFactor[\"$pro\"] = 1.0  # $pro"
  end
  # cooperativity number in transfer function
  buffer *= "\n\n\t# cooperativity number in transfer function"
  buffer *= "\n\tcooperativity = 1"
  #buffer *= "\n\n\tbackground_mRNA_synthesis_rate_vector = 0.01*ones($length_TXTL)"
  # load txtl constants buffer
  buffer *= "\n\n"
  if host_type == "bacteria"
    buffer *= jl_include_constants_from_literature(
              joinpath(Base.@__DIR__, "txtl_constants_ecoli.jl"),"\n\t")
  else
    buffer *= jl_include_constants_from_literature(
              joinpath(Base.@__DIR__, "txtl_constants_hl60.jl"), "\n\t")
  end

  #---------------------------------
  # put all stuff in a Dictionary
  buffer *= "\n\n\t####################################"
  buffer *= "\n\t# put all stuff in a Dictionary"
  buffer *= "\n\tdataDictionary = Dict{String, Any}()"
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
  buffer *= "\nend"

  return buffer
end

function jl_build_simulation_buffer(NoExtracellularSpecies::Int64)
  # buffer = build_copyright_header_buffer()
  buffer = "\n# set up ODE, get the derivatives
function Balances(t,y,dataDictionary)
  # correct for negatives
  idx_small = findall(y.<0)
  y[idx_small] .= 0.0

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
  (rnx_rate, transcription_rate, translation_rate) = calculate_kinetics(y, dataDictionary)

  # calculate the derivatives
  dydt = zeros(length(all_species_dict))

  # signaling reaction network
  X2 = zeros(length(rnx_species)-$(NoExtracellularSpecies))  # species inside the cell
  for (id, st) in enumerate(rnx_species[$(NoExtracellularSpecies)+1:end])  # get X2 value from y
    X2[id] = y[all_species_dict[st]]
  end
  dX1 = stoichiometry[1:$(NoExtracellularSpecies), :]*rnx_rate  # extracellular metabolites
  dX2 = stoichiometry[$(NoExtracellularSpecies)+1:end, :]*rnx_rate - specific_growth_rate*X2  # INTRA
  dX = [dX1; dX2]
  for (id, st) in enumerate(rnx_species)  # return derivatives back to dydt
    dydt[all_species_dict[st]] = dX[id]
  end
  dydt[1] = specific_growth_rate*y[1]  # for Biomass

  # TXTL network
  for i = 1:length(mRNA_species)
  # seems easy to add additional terms to account for the signaling effect
    mRNA_id = all_species_dict[mRNA_species[i]]
    p_id = all_species_dict[protein_species[i]]
    dydt[mRNA_id] = transcription_rate[mRNA_species[i]] - (specific_growth_rate + transcription_correction[mRNA_species[i]]*degradation_constant_mRNA) * y[mRNA_id]
    dydt[p_id] = translation_rate[protein_species[i]] - (specific_growth_rate + translation_correction[protein_species[i]]*degradation_constant_protein) * y[p_id]
  end

  return dydt
end"
  return buffer
end


function jl_build_solveODEBalances_buffer(all_species_array::Array,
  all_species_dict::Dict{String, Int}, mRNA_species::Array, protein_species::Array)

  buffer = "include(\"./Include.jl\")" *
    "\n\nfunction SolveBalances(TStart, TStop, TStep, dataDictionary)" *
    "\n\tt = collect(TStart:TStep:TStop)" *
    "\n\tall_species_reversed_dict = dataDictionary[\"all_species_reversed_dict\"]"
  # initial concentration
  buffer *= "\n\t# initial species concentration"
  buffer *= "\n\ty0 = dataDictionary[\"initial_condition\"]"
  # run simulation
  buffer *= "\n\t# call ODE solver"
  buffer *= "\n\tf(t, y) = Balances(t, y, dataDictionary)" *
    "\n\tt, y = ode23s(f, y0, t; points=:specified)"
  # data transfer
  buffer *= "\n\n\t# data transfer" *
    "\n\trow = length(t)" *
    "\n\tcol = length(y0)" *
    "\n\tY = zeros(row, col)" *
    "\n\tforeach(x->(Y[x, :] = y[x]), collect(1:row))"
  # plotting
  buffer *= "\n\t# plot results" *
    "\n\tplt.figure(\"simulation \$TStart to \$TStop\")" *
    "\n\tsubplt_col = ceil(Int, col/4)" *
    "\n\tfor i = 1:col" *
    "\n\t\tplt.subplot(4, subplt_col, i)" *
    "\n\t\tplt.plot(t, Y[:,i])" *
    "\n\t\tplt.title(all_species_reversed_dict[i])" *
    "\n\tend" *
    "\n\tplt.show()"
  buffer *= "\n\n\treturn (t, Y)" *
    "\nend"

  # add a test
  buffer *= "\n\n\n"
  buffer *= "# # uncomment to test"
  buffer *= "\n# dataDictionary = generate_model_parameters_dictionary()"
  buffer *= "\n# SolveBalances(0, 10, 0.1, dataDictionary)"

  return buffer
end

# function build_copyright_header_buffer()
#   current_year = string(Dates.year(now()))
#   # Get comment data from
#   buffer = ""
#   buffer *=
# "# ----------------------------------------------------------------------------------- #
# # Copyright (c) $(current_year) Varnerlab
# # Robert Frederick Smith School of Chemical and Biomolecular Engineering
# # Cornell University, Ithaca NY 14850
# #
# # Permission is hereby granted, free of charge, to any person obtaining a copy
# # of this software and associated documentation files (the \"Software\"), to deal
# # in the Software without restriction, including without limitation the rights
# # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# # copies of the Software, and to permit persons to whom the Software is
# # furnished to do so, subject to the following conditions:
# #
# # The above copyright notice and this permission notice shall be included in
# # all copies or substantial portions of the Software.
# #
# # THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# # THE SOFTWARE.
# # ----------------------------------------------------------------------------------- #"
#
#   return buffer
# end


# FBA data dictionary generation
function jl_generate_FBA_data_dictionary(all_rnx_list::Array,
  rnx_species_array::Array, extra_species_num::Int)
  secrete_id = Set{Int64}()  # For initialize the coefficient array
  # data dictionary
  buffer = "function DataDictionary(time_start,time_stop,time_step)"
  # stoichiometry
  buffer *= "\n\n\t# load stoichiometry"
  buffer *= "\n\tstoichiometric_matrix = readdlm(\"./stoichiometry.dat\")"
  # Setup default flux bounds array
  buffer *= "\n\t# Setup default flux bounds array"
  buffer *= "\n\tdefault_bounds_array = ["
  for (id, rnx) in enumerate(all_rnx_list)
    buffer *= "\n\t\t0 100.0; # $(id) $(rnx.rnxName)"
  end
  buffer *= "\n\t]"
  # Setup default species bounds array
  buffer *= "\n\n\t# Setup default species bounds array"
  buffer *= "\n\tspecies_bounds_array = ["
  for (id, sp) in enumerate(rnx_species_array)
    if id > extra_species_num
      buffer *= "\n\t\t 0.0 0.0; # $(id) $(sp)"
    else
      buffer *= "\n\t\t-1.0 1.0; # $(id) $(sp)"
    end
  end
  buffer *= "\n\t]"
  # Setup the objective coefficient array
  buffer *= "\n\n\t# Setup the objective coefficient array"
  buffer *= "\n\tobjective_coefficient_array = ["
  for (id, rnx) in enumerate(all_rnx_list)
    buffer *= "\n\t\t0.0; # $(id) $(rnx.rnxName)"
    if rnx.rnxType == "secrete"
        push!(secrete_id, id)
    end
  end
  buffer *= "\n\t]"
  # Min/Max flag - default is minimum -
  buffer *= "\n\n\t# Min/Max flag - default is minimum -"
  buffer *= "\n\tis_minimum_flag = true"
  # List of reation strings - used to write flux report
  buffer *= "\n\n\t# List of reation strings - used to write flux report"
  buffer *= "\n\tlist_of_reaction_strings = ["
  for rnx in all_rnx_list
    buffer *= "\n\t\t\"$(rnx.rnxName)\""
  end
  buffer *= "\n\t]"
  # List of metabolite strings - used to write flux report
  buffer *= "\n\n\t# List of metabolite strings - used to write flux report"
  buffer *= "\n\tlist_of_metabolite_symbols = ["
  for sp in rnx_species_array
    buffer *= "\n\t\t\"$(sp)\""
  end
  buffer *= "\n\t]"
  # return block -
  buffer *= "\n\n"
  buffer *= "\t# =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
  buffer *= "\tdata_dictionary = Dict{AbstractString,Any}()\n"
  buffer *= "\tdata_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix\n"
  buffer *= "\tdata_dictionary[\"objective_coefficient_array\"] = objective_coefficient_array\n"
  buffer *= "\tdata_dictionary[\"default_flux_bounds_array\"] = default_bounds_array;\n"
  buffer *= "\tdata_dictionary[\"species_bounds_array\"] = species_bounds_array\n"
  buffer *= "\tdata_dictionary[\"list_of_reaction_strings\"] = list_of_reaction_strings\n"
  buffer *= "\tdata_dictionary[\"list_of_metabolite_symbols\"] = list_of_metabolite_symbols\n"
  buffer *= "\tdata_dictionary[\"is_minimum_flag\"] = is_minimum_flag\n"
  buffer *= "\t# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
  buffer *= "\treturn data_dictionary\n"
  buffer *= "end\n"

  # maximize_product_dictionary
  buffer *= "\n\n\n# maximize_product_dictionary"
  buffer *= "\nfunction maximize_product_dictionary(time_start,time_stop,time_step)"
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
  buffer *= "\nend"

  return buffer
end


# Kinetic model Include.jl
function jl_build_kinetics_include_buffer()
 buffer =
"include(\"DataDictionary.jl\")
include(\"Kinetics.jl\")
include(\"Balances.jl\")

import PyPlot
using ODE
using DelimitedFiles

const plt = PyPlot"
  return buffer
end

# FBA/FVA model Include.jl
function jl_build_FBA_include_buffer()
 buffer =
"include(\"DataDictionary.jl\")
include(\"FluxDriver.jl\")
include(\"FVA.jl\")

using GLPK
using DelimitedFiles"
  return buffer
end



# FBA/FVA Solve.jl
function jl_build_FBA_solve_buffer()
  buffer =
"include(\"Include.jl\")

# load the data dictionary -
data_dictionary = maximize_product_dictionary(0,10,1)

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
(calculated_flux_array, dual_value_array) = calculate_flux_variabilty(data_dictionary)

# FBA
println(\"objective: \$objective_value\")
println(\"flux array: \", flux_array)
println(\"uptake array: \", uptake_array)

# FVA
println(\"FVA results: \")
for i = 1:size(calculated_flux_array)[1]
  println(\"    \", calculated_flux_array[i,:])
end"
  return buffer
end


# rFBA model Include.jl
function jl_build_rFBA_include_buffer()
 buffer =
"include(\"DataDictionary.jl\")
include(\"FluxDriver.jl\")
include(\"FVA.jl\")
include(\"rRules.jl\")

using GLPK
using DelimitedFiles"
  return buffer
end


# rFBA rRules.jl
function jl_build_rRules_buffer()
  buffer =
"function rRules(bounds, args...)
  # bounds: lower/upper bounds on flux
  # args: variable number of arguments, optional

  # initialize regulatory factors
  rfs = ones(size(bounds))
  # fill in regulatory rules below

  # update bounds
  new_bounds = rfs.*bounds

  return new_bounds
end"
  return buffer
end


# rFBA Solve.jl
function jl_build_rFBA_solve_buffer()
  buffer =
"include(\"Include.jl\")

# load the data dictionary -
data_dictionary = maximize_product_dictionary(0,10,1)

# apply regulatory rules
new_bounds = rRules(data_dictionary[\"default_flux_bounds_array\"])
data_dictionary[\"default_flux_bounds_array\"] = new_bounds

# solve the lp problem -
(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)
(calculated_flux_array, dual_value_array) = calculate_flux_variabilty(data_dictionary)

# FBA
println(\"objective: \$objective_value\")
println(\"flux array: \", flux_array)
println(\"uptake array: \", uptake_array)

# FVA
println(\"FVA results: \")
for i = 1:size(calculated_flux_array)[1]
  println(\"    \", calculated_flux_array[i,:])
end"
  return buffer
end
