
# set up ODE, get the derivatives
function Balances(t,y,dataDictionary)
  # correct for negatives
  idx_small = findall(y.<0)
  y[idx_small] .= 0.0

  # load data from data dictionary
  rnx_species = dataDictionary["rnx_species_array"]
  mRNA_species = dataDictionary["mRNA_species_array"]
  protein_species = dataDictionary["protein_species_array"]
  stoichiometry = dataDictionary["stoichiometric_matrix"]
  all_species_dict = dataDictionary["all_species_dict"]  # String --> number
  specific_growth_rate = dataDictionary["specificGrowthRate"]
  transcription_correction = dataDictionary["transcriptionSpecificCorrectionFactorDict"]
  translation_correction = dataDictionary["translationSpecificCorrectionFactor"]
  degradation_constant_mRNA = dataDictionary["degradationConstantmRNA"]
  degradation_constant_protein = dataDictionary["degradationConstantProtein"]

  # kinetics
  (rnx_rate, transcription_rate, translation_rate) = calculate_kinetics(y, dataDictionary)

  # calculate the derivatives
  dydt = zeros(length(all_species_dict))

  # signaling reaction network
  X2 = zeros(length(rnx_species)-3)  # species inside the cell
  for (id, st) in enumerate(rnx_species[3+1:end])  # get X2 value from y
    X2[id] = y[all_species_dict[st]]
  end
  dX1 = stoichiometry[1:3, :]*rnx_rate  # extracellular metabolites
  dX2 = stoichiometry[3+1:end, :]*rnx_rate - specific_growth_rate*X2  # INTRA
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
end