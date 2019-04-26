function calculate_kinetics(X, data_dictionary)

	# load all species dictionary
	all_species_dict = data_dictionary["all_species_dict"]  # String-->Int

	# generate kinetics for signaling, assume saturation kinetics
	# load Monod affinity constants
	MonodK = data_dictionary["MonodAffinityConstantDict"]  # String-->Float
	# load reaction kinetic constants
	kcat = data_dictionary["kcat_signaling"]  # in order of rnx

	# write reaction rate equations
	rnx_rate_vector = zeros(7)
	rnx_rate_vector[1] = kcat[1]*X[12]*X[2]/(MonodK["MonodK~rnx1~m_A_e"]+X[2])
	rnx_rate_vector[2] = kcat[2]*X[13]*X[6]/(MonodK["MonodK~rnx2~m_B_c"]+X[6])
	rnx_rate_vector[3] = kcat[3]*X[14]*X[7]/(MonodK["MonodK~rnx3~m_C_c"]+X[7])
	rnx_rate_vector[4] = kcat[4]*X[8]*X[5]/(MonodK["MonodK~rnx4~m_A_c"]+X[5])
	rnx_rate_vector[5] = kcat[5]*X[11]*X[6]/(MonodK["MonodK~rnx5~m_B_c"]+X[6])
	rnx_rate_vector[6] = kcat[6]*X[9]*X[5]/(MonodK["MonodK~rnx6~m_A_c"]+X[5])
	rnx_rate_vector[7] = kcat[7]*X[10]*X[7]/(MonodK["MonodK~rnx7~m_C_c"]+X[7])

	# generate kinetics & control terms for TXTL
	# load data from data dictionary
	W_value_dict = data_dictionary["W_value_dict"]  #String-->Float
	coop = data_dictionary["cooperativity"]
	disassociation_constant_dict = data_dictionary["transferFunctionDisassociationConstantDict"]
	background_control_term = data_dictionary["backgroundGeneExpressionControlTermDict"]
	kcatTX = data_dictionary["kcatTranscription"]
	RNAPconc = data_dictionary["RNAPConcentration"]
	geneConc = data_dictionary["avgGeneConcentration"]
	saturationTX = data_dictionary["transcriptionSaturationConstant"]
	kcatTL = data_dictionary["kcatTranslation"]
	RIBOconc = data_dictionary["RIBOConcentration"]
	saturationTL = data_dictionary["translationSaturationConstant"]

	# write control terms, transcription rate equations, and translation rate equations
	TX_control_term = Dict{String, Float64}()  # control term
	TX_rate_vector = Dict{String, Float64}()  # transcription rate
	TL_rate_vector = Dict{String, Float64}()  # translation rate
	# mRNA_E4_c and p_E4_c
	mRNA_E4_c_m_B_c = W_value_dict["W~m_B_c~mRNA_E4_c"] * X[6]^coop/((disassociation_constant_dict["KD~m_B_c"])^coop+X[6]^coop)
	up_action_on_mRNA_E4_c = mRNA_E4_c_m_B_c 
	TX_control_term["mRNA_E4_c"] = (background_control_term["mRNA_E4_c"] + up_action_on_mRNA_E4_c)/(1 + background_control_term["mRNA_E4_c"] + up_action_on_mRNA_E4_c)
	TX_rate_vector["mRNA_E4_c"] = TX_control_term["mRNA_E4_c"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)
	TL_rate_vector["p_E4_c"] = kcatTL*RIBOconc*X[18]/(saturationTL+X[18])
	# mRNA_E2_c and p_E2_c
	mRNA_E2_c_m_A_c = W_value_dict["W~m_A_c~mRNA_E2_c"] * X[5]^coop/((disassociation_constant_dict["KD~m_A_c"])^coop+X[5]^coop)
	up_action_on_mRNA_E2_c = mRNA_E2_c_m_A_c 
	TX_control_term["mRNA_E2_c"] = (background_control_term["mRNA_E2_c"] + up_action_on_mRNA_E2_c)/(1 + background_control_term["mRNA_E2_c"] + up_action_on_mRNA_E2_c)
	TX_rate_vector["mRNA_E2_c"] = TX_control_term["mRNA_E2_c"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)
	TL_rate_vector["p_E2_c"] = kcatTL*RIBOconc*X[16]/(saturationTL+X[16])
	# mRNA_E3_c and p_E3_c
	mRNA_E3_c_m_C_c = W_value_dict["W~m_C_c~mRNA_E3_c"] * X[7]^coop/((disassociation_constant_dict["KD~m_C_c"])^coop+X[7]^coop)
	up_action_on_mRNA_E3_c = mRNA_E3_c_m_C_c 
	TX_control_term["mRNA_E3_c"] = (background_control_term["mRNA_E3_c"] + up_action_on_mRNA_E3_c)/(1 + background_control_term["mRNA_E3_c"] + up_action_on_mRNA_E3_c)
	TX_rate_vector["mRNA_E3_c"] = TX_control_term["mRNA_E3_c"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)
	TL_rate_vector["p_E3_c"] = kcatTL*RIBOconc*X[17]/(saturationTL+X[17])
	# mRNA_TA_c and p_TA_c
	mRNA_TA_c_m_A_e = W_value_dict["W~m_A_e~mRNA_TA_c"] * X[2]^coop/((disassociation_constant_dict["KD~m_A_e"])^coop+X[2]^coop)
	up_action_on_mRNA_TA_c = mRNA_TA_c_m_A_e 
	TX_control_term["mRNA_TA_c"] = (background_control_term["mRNA_TA_c"] + up_action_on_mRNA_TA_c)/(1 + background_control_term["mRNA_TA_c"] + up_action_on_mRNA_TA_c)
	TX_rate_vector["mRNA_TA_c"] = TX_control_term["mRNA_TA_c"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)
	TL_rate_vector["p_TA_c"] = kcatTL*RIBOconc*X[19]/(saturationTL+X[19])
	# mRNA_TB_c and p_TB_c
	mRNA_TB_c_m_B_c = W_value_dict["W~m_B_c~mRNA_TB_c"] * X[6]^coop/((disassociation_constant_dict["KD~m_B_c"])^coop+X[6]^coop)
	up_action_on_mRNA_TB_c = mRNA_TB_c_m_B_c 
	TX_control_term["mRNA_TB_c"] = (background_control_term["mRNA_TB_c"] + up_action_on_mRNA_TB_c)/(1 + background_control_term["mRNA_TB_c"] + up_action_on_mRNA_TB_c)
	TX_rate_vector["mRNA_TB_c"] = TX_control_term["mRNA_TB_c"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)
	TL_rate_vector["p_TB_c"] = kcatTL*RIBOconc*X[20]/(saturationTL+X[20])
	# mRNA_E1_c and p_E1_c
	mRNA_E1_c_m_A_c = W_value_dict["W~m_A_c~mRNA_E1_c"] * X[5]^coop/((disassociation_constant_dict["KD~m_A_c"])^coop+X[5]^coop)
	up_action_on_mRNA_E1_c = mRNA_E1_c_m_A_c 
	TX_control_term["mRNA_E1_c"] = (background_control_term["mRNA_E1_c"] + up_action_on_mRNA_E1_c)/(1 + background_control_term["mRNA_E1_c"] + up_action_on_mRNA_E1_c)
	TX_rate_vector["mRNA_E1_c"] = TX_control_term["mRNA_E1_c"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)
	TL_rate_vector["p_E1_c"] = kcatTL*RIBOconc*X[15]/(saturationTL+X[15])
	# mRNA_TC_c and p_TC_c
	mRNA_TC_c_m_C_c = W_value_dict["W~m_C_c~mRNA_TC_c"] * X[7]^coop/((disassociation_constant_dict["KD~m_C_c"])^coop+X[7]^coop)
	up_action_on_mRNA_TC_c = mRNA_TC_c_m_C_c 
	TX_control_term["mRNA_TC_c"] = (background_control_term["mRNA_TC_c"] + up_action_on_mRNA_TC_c)/(1 + background_control_term["mRNA_TC_c"] + up_action_on_mRNA_TC_c)
	TX_rate_vector["mRNA_TC_c"] = TX_control_term["mRNA_TC_c"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)
	TL_rate_vector["p_TC_c"] = kcatTL*RIBOconc*X[21]/(saturationTL+X[21])

	return (rnx_rate_vector, TX_rate_vector, TL_rate_vector)
end