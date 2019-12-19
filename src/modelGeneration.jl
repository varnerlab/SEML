# model generation
mutable struct reactionForm
  reactants::Array
  products::Array
  catalysts::Array
  rnxName::AbstractString
  rnxType::AbstractString
  function reactionForm()
    new()
  end
end
mutable struct txtlForm
  activationProtein::Array
  inhibitionProtein::Array
  function txtlForm()
    new([], [])
  end
end

#=
F: collect reactions in normal form, also get rid of "[]"
I: Array of tuples, tuple = (reaction type, [reactants, products, catalysts]).
 [reactants, products, catalysts] are all in List form
O: Array of reactionForm, Array of txtlForm,
  Set of rnx species, Set of mRNA species, Set of protein species from mRNA
BA
=#
function preparing_rnxList_txtlDict(preIRtuples::Array, sys2user::Dict, user2sys::Dict)
  # preIRtuples = array of tuple ((verb, line #), [reactant, product, catalyst])
  rnxList = Array{reactionForm,1}()
  rnx_species_set = Set{String}()

  txtlDict = Dict{AbstractString, txtlForm}()
  txtl_mRNA_set = Set{String}()
  txtl_protein_set = Set{String}()
  # generating rnx but also remember to add species to set
  for tuple in preIRtuples
    if tuple[1][1] == "react"  # get rid of []
      tmpRF = reactionForm()
      tmpRF.rnxType = tuple[1][1]
      name_r = ""
      name_c = ""
      name_p = ""
      # println(tuple[2][1])
      if tuple[2][1][1].oriBioName == "[]"
        tmpRF.products = tuple[2][2]
        name_p = join(map(s->string(s.coeff)*"*"*string(s.oriBioName), tmpRF.products), "+")
        foreach(s->push!(rnx_species_set, s.oriBioName), tmpRF.products)
      elseif tuple[2][2][1].oriBioName == "[]"
        tmpRF.reactants = tuple[2][1]
        name_r = join(map(s->string(s.coeff)*"*"*string(s.oriBioName), tmpRF.reactants), "+")
        foreach(s->push!(rnx_species_set, s.oriBioName), tmpRF.reactants)
      else
        println("WARNING: ERROR in \'react\', \'[]\' is not processed properly")
      end
      tmpRF.rnxName = name_r *"<"* tmpRF.rnxType *":"* name_c *">"* name_p
      push!(rnxList, tmpRF)
    elseif (tuple[1][1] == "uptake" || tuple[1][1] == "secrete" || tuple[1][1] == "bind" ||
      tuple[1][1] == "unbind" || tuple[1][1] == "phosphorylate" ||
      tuple[1][1] == "dephosphorylate" || tuple[1][1] == "catalyze")
      tmpRF = reactionForm()
      tmpRF.rnxType = tuple[1][1]
      name_r = ""
      name_c = ""
      name_p = ""
      tmpRF.reactants = tuple[2][1]
      name_r = join(map(s->string(s.coeff)*"*"*string(s.oriBioName), tmpRF.reactants), "+")
      foreach(s->push!(rnx_species_set, s.oriBioName), tmpRF.reactants)
      tmpRF.products = tuple[2][2]
      name_p = join(map(s->string(s.coeff)*"*"*string(s.oriBioName), tmpRF.products), "+")
      foreach(s->push!(rnx_species_set, s.oriBioName), tmpRF.products)
      if isassigned(tuple[2], 3) && tuple[2][3][1].oriBioName != "[]" # has catalysts
        tmpRF.catalysts = tuple[2][3]
        name_c = join(map(s->string(s.coeff)*"*"*string(s.oriBioName), tmpRF.catalysts), "+")
        foreach(s->push!(rnx_species_set, s.oriBioName), tmpRF.catalysts)
      end
      tmpRF.rnxName = name_r *"<"* tmpRF.rnxType *":"* name_c *">"* name_p
      push!(rnxList, tmpRF)

    elseif tuple[1][1] == "induce"
      process_txtl_relation(txtlDict, tuple[2], sys2user["GENE"], sys2user["MRNA"], 1, txtlForm)
    elseif tuple[1][1] == "activate"
      process_txtl_relation(txtlDict, tuple[2], sys2user["PROTEIN"], sys2user["MRNA"], 1, txtlForm)
    elseif tuple[1][1] == "repress"
      process_txtl_relation(txtlDict, tuple[2], sys2user["GENE"], sys2user["MRNA"], 2, txtlForm)
    elseif tuple[1][1] == "inhibit"
      process_txtl_relation(txtlDict, tuple[2], sys2user["PROTEIN"], sys2user["MRNA"], 2, txtlForm)
    end
  end

  txtl_mRNA_set = Set(keys(txtlDict))
  # pre_txtl_protein_set = [replace(s, sys2user["MRNA"] => sys2user["PROTEIN"]) for s in txtl_mRNA_set]
  txtl_protein_set = Set([replace(s, sys2user["MRNA"] => sys2user["PROTEIN"]) for s in txtl_mRNA_set])


  return (rnxList, txtlDict, rnx_species_set, txtl_mRNA_set, txtl_protein_set)
end

#=
F: add gene regulators;
I:
O: modified txtl;
BA
=#
function process_txtl_relation(txtl::Dict, tuple::Array, oldTag::AbstractString,
  newTag::AbstractString, fieldNo::Int, dataform::Type)
  targetMRNA = replace(tuple[2][1].oriBioName, oldTag => newTag)
  if !haskey(txtl, targetMRNA)  # no exist -> create it
    txtl[targetMRNA] = dataform()
  end
  if fieldNo == 1
    push!(txtl[targetMRNA].activationProtein, tuple[1])
  else
    push!(txtl[targetMRNA].inhibitionProtein, tuple[1])
  end
end

#=
F: sort species symbols for convenience;
I:
O: ;
B: how to support user-defined delimiters for type and compartment? e.g. _exc,
  -exc, or ~ext? and how to guarantee that compartment tag is always at the end
  while type tag is the first?
A:
=#
function sorting_species_list(unsorted_species_set::Set, mRNA_set::Set,
  mprotein_set::Set, sys2user::Dict)
  extracellular_metabolite_array = Array{String,1}()
  metabolite_array = Array{String,1}()
  protein_array = Array{String,1}()
  mRNA_array = Array{String,1}()
  the_else_array = Array{String,1}()
  # determing delimiter  FIXME
  if occursin("_", sys2user["_exc"])
    delimiter_com = "_"
  elseif occursin("-", sys2user["_exc"])
    delimiter_com = "-"
  elseif occursin("~", sys2user["_exc"])
    delimiter_com = "~"
  else
    delimiter_com = "_"
    println("WARNING: weired compartment tag $(sys2user["_exc"])")
  end
  if occursin("_", sys2user["PROTEIN"])
    delimiter_type = "_"
  elseif occursin("-", sys2user["PROTEIN"])
    delimiter_type = "-"
  elseif occursin("~", sys2user["PROTEIN"])
    delimiter_type = "~"
  else
    delimiter_type = "_"
    println("WARNING: weired type tag $(sys2user["PROTEIN"])")
  end
  # sort by type
  for st in unsorted_species_set
    st_pieces = split(st, delimiter_com)
    st_compartment = delimiter_com*st_pieces[end]
    st_p2 = split(st, delimiter_type)
    st_type = st_p2[1]*delimiter_type
    if st_compartment == sys2user["_exc"]
      push!(extracellular_metabolite_array, st)
    elseif st_type == sys2user["METABOLITE"]
      push!(metabolite_array, st)
    elseif st_type == sys2user["PROTEIN"]
      push!(protein_array, st)
    elseif st_type == sys2user["MRNA"]
      push!(mRNA_array, st)
    else  # weired symbol?
      push!(the_else_array, st)
      println("WARNING: weired bio-symbol: $(st)")
    end
  end
  # processing results
  sorted_species_array = [sort(extracellular_metabolite_array); sort(metabolite_array);
    sort(protein_array); sort(mRNA_array); sort(the_else_array)]
  species2index_dict = Dict{String, Int64}()  # create species dict: string --> int64
  index2species_dict = Dict{Int64, String}()  # create species dict: int64 --> string
  for (index, st) in enumerate(sorted_species_array)
    species2index_dict[st] = index
    index2species_dict[index] = st
  end
  # merge results
  sorted_all_species_array = ["BIOMASS"; sort(extracellular_metabolite_array);
    sort(metabolite_array); sort(collect(union(Set(protein_array), mprotein_set)));
    sort(collect(union(Set(mRNA_array), mRNA_set))); sort(the_else_array)]
  all_species2index_dict = Dict{String, Int64}()  # create species dict: string --> int64
  all_index2species_dict = Dict{Int64, String}()  # create species dict: int64 --> string
  for (index, st) in enumerate(sorted_all_species_array)
    all_species2index_dict[st] = index
    all_index2species_dict[index] = st
  end

  return (sorted_species_array, species2index_dict, index2species_dict,
    sorted_all_species_array, all_species2index_dict, all_index2species_dict,
    length(extracellular_metabolite_array))
end

#=
F: stoichiometry;
I:
O: ;
B:
A:
=#
function build_stoichiometric_matrix_buffer(all_reaction_list::Array, rnx_species_dict::Dict{String, Int64})
  stoichiometric_matrix = zeros(length(rnx_species_dict), length(all_reaction_list)) # generate matrix
  for (index, rnx) in enumerate(all_reaction_list)  # go thru each rnx
    if isdefined(rnx, :reactants)
      foreach(s->(stoichiometric_matrix[rnx_species_dict[s.oriBioName], index] = -s.coeff), rnx.reactants)
    end
    if isdefined(rnx,:products)
      foreach(s->(stoichiometric_matrix[rnx_species_dict[s.oriBioName], index] = s.coeff), rnx.products)
    end
  end
  (row, col) = size(stoichiometric_matrix)
  stoichiometric_buffer = ""  # convert into string
  for i = 1:row
    for j = 1:col
      stoichiometric_buffer *= "$(stoichiometric_matrix[i,j]) "
    end
    stoichiometric_buffer = chop(stoichiometric_buffer)
    stoichiometric_buffer *= "\n"
  end
  return stoichiometric_buffer, stoichiometric_matrix
end
