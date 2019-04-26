# washout simulation
include("Include.jl")
include("SolveBalances.jl")

# set up simulation time for each phase
time_start = 0
PI_dura = 7
PII_dura = 10
PIII_dura = 8
time_step = 0.1

# phase I: run to steady state
DataDict = generate_model_parameters_dictionary()
TIEnd = time_start + PI_dura
(TP1, XP1) = SolveBalances(time_start, TIEnd, time_step, DataDict)
# phase II: add substrate
IniCond = XP1[end, :]
IniCond[1] = 0.001
IniCond[2] = 10
DataDict["initial_condition"] = IniCond
TIIEnd = TIEnd + PII_dura
(TP2, XP2) = SolveBalances(TIEnd, TIIEnd, time_step, DataDict)
# phase III: remove substrate
IniCond = XP2[end, :]
IniCond[2] = 0
DataDict["initial_condition"] = IniCond
TIIIEnd = TIIEnd + PIII_dura
(TP3, XP3) = SolveBalances(TIIEnd, TIIIEnd, time_step, DataDict)

# pack the three phases together
T = [TP1; TP2; TP3]
X = [XP1 ; XP2 ; XP3];
# plot results
all_species_reversed_dict = DataDict["all_species_reversed_dict"]

# metG = 3:7
# proG = 8:14
# mRNAG = 15:21
# group = ["metabolite", "protein", "mRNA"]
# sDt = Dict("metabolite" => metG, "protein" => proG, "mRNA" => mRNAG)
# for st in group
#   plt.figure(st)
#   # plt.plot(T[40:end]-4, X[40:end,sDt[st]], label=[all_species_reversed_dict[i] for i in sDt[st]])
#   for i in sDt[st]
#     plt.plot(T[40:end]-4, X[40:end,i], label=all_species_reversed_dict[i])
#   end
#   plt.title(st)
#   plt.legend()
#   plt.show()
# end


function plot_by_groups(groups, labels)
  for (key,val) in groups
    plt.figure(key)
    for i in val
      plt.plot(T[40:end] .- 4, X[40:end,i], label=all_species_reversed_dict[i])
    end
    plt.axvspan(0, 3, color="k", alpha=0.1)
    plt.axvspan(3,13, color="g", alpha=0.1)
    plt.axvspan(13,21, color="b", alpha=0.1)
    plt.title(key)
    plt.xlim(0, 21)
    plt.ylim(bottom=0)
    plt.xlabel("Time (h)")
    plt.ylabel("Concentration ($(labels[key]))")
    plt.legend()
    # plt.show()
  end
end

sDt = Dict("metabolite" => 3:7, "protein" => 8:14, "mRNA" => 15:21)
ylab = Dict("metabolite" => "mM", "protein" => "nM", "mRNA" => "nM")
plot_by_groups(sDt, ylab)
plt.show()
# sDt2 = Dict("TA" => [19,12,5], "A" => [5,15,16,8,9,6,7],
#             "C" => [7,17,21,10,14,4,6], "B" => [6,18,20,11,13,6])
# plot_by_groups(sDt2)

# writedlm("simulation.dat", [T X])

# column; species
# 1  BIOMASS
# 2  met_A_e
# 3  met_B_e
# 4  met_C_e
# 5  met_A_c
# 6  met_B_c
# 7  met_C_c
# 8  p_E1_c
# 9  p_E2_c
# 10  p_E3_c
# 11  p_E4_c
# 12  p_TA_c
# 13  p_TB_c
# 14  p_TC_c
# 15  mRNA_E1_c
# 16  mRNA_E2_c
# 17  mRNA_E3_c
# 18  mRNA_E4_c
# 19  mRNA_TA_c
# 20  mRNA_TB_c
# 21  mRNA_TC_c
