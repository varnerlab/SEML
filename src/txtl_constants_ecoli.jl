# ------------------------------------------------------------------------------------------#
# constants (from bionumbers)       units
# ------------------------------------------------------------------------------------------#
cell_diameter = 1.1                 # mum
number_of_rnapII = 4600            	# copies/cells
number_of_ribosome = 50000         	# copies/cells
mRNA_half_life_TF = 0.083           # hrs
protein_half_life = 70              # hrs
doubling_time_cell = 0.33           # hrs
max_translation_rate = 16.5         # aa/sec
max_transcription_rate = 60.0       # nt/sec
average_transcript_length = 1200   	# nt
average_protein_length = 400       	# aa
fraction_nucleus = 0.0             	# dimensionless
av_number = 6.02e23                 # number/mol
avg_gene_number = 2                 # number of copies of a gene
polysome_number = 4									# number of ribsomoses per transcript
# ------------------------------------------------------------------------------------------#
#
# ------------------------------------------------------------------------------------------#
# Calculate constants using bionumber values
# ------------------------------------------------------------------------------------------#
# Calculate the volume (convert to L)
V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)     # volume

# Calculate the rnapII_concentration and ribosome_concentration
rnapII_concentration = number_of_rnapII*(1/av_number)*(1/V)*1e9                   # nM
ribosome_concentration = number_of_ribosome*(1/av_number)*(1/V)*1e9               # nM

# degrdation rate constants -
degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                       # hr^-1
degradation_constant_protein = -(1/protein_half_life)*log(0.5)                    # hr^-1

# kcats for transcription and translation -
kcat_transcription = max_transcription_rate*(3600/average_transcript_length)      # hr^-1
kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)             # hr^-1

# Maximum specific growth rate -
maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                      # hr^-1

# What is the average gene concentration -
avg_gene_concentration = avg_gene_number*(1/av_number)*(1/V)*1e9                  # nM

# How fast do my cells die?
death_rate_constant = 0.05*maximum_specific_growth_rate                            # hr^-1

# Saturation constants for translation and trascription -
saturation_transcription = 4600*(1/av_number)*(1/V)*1e9                           # nM
saturation_translation = 150000*(1/av_number)*(1/V)*1e9                           # nM
# -------------------------------------------------------------------------------------------#
