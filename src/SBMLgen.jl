using PyCall

function generate_SBML_file(sys2userDict::Dict, rnx_set::Set, rnx::Array)
  sbml = pyimport("libsbml")
  sbmlns = sbml.SBMLNamespaces(3,2)
  #  create the document
  document = sbml.SBMLDocument(sbmlns)
  #  create the Model
  model= document.createModel()

  # -----helper functions start------
  #  create the Compartment
  function create_compartment(m, ID)
    compartment = m.createCompartment()
    compartment.setId(ID)
    # compartment.setConstant(true)
    # compartment.setSize(1)
  end
  # create the Species
  function create_species(m, ID, cpt)
    species = m.createSpecies()
    species.setId(ID)
    species.setCompartment(cpt)
    # species.setBoundaryCondition(false)
    # species.setConstant(false)
  end
  # create reactant
  function create_reactant(rn, ID, coe)
    reactant = rn.createReactant()
    reactant.setSpecies(ID)
    reactant.setStoichiometry(coe)
  end
  # create product
  function create_product(rn, ID, coe)
    product = rn.createProduct()
    product.setSpecies(ID)
    product.setStoichiometry(coe)
  end
  # create reactions
  function create_reaction(m, rxn)
    reaction = m.createReaction()
    Id = replace(replace(rxn.rnxName, r"[0-9]|\.|\*" => ""), r"<|:>|>|:" => "__")
    reaction.setId(Id)
    reaction.setReversible(false)
    for rt in rxn.reactants
      create_reactant(reaction, rt.oriBioName, rt.coeff)
    end
    for pt in rxn.products
      create_product(reaction, pt.oriBioName, pt.coeff)
    end
  end
  # -----helper functions end------

  # compartments
  if haskey(sys2userDict, "_exc")
    create_compartment(model, sys2userDict["_exc"])
  end
  if haskey(sys2userDict, "_c")
    create_compartment(model, sys2userDict["_c"])
  end
  # species
  for sps in rnx_set
    create_species(model, sps, "_" * split(sps, "_")[end])
  end
  # reactions
  for rxn in rnx
    create_reaction(model, rxn)
  end

  return sbml.writeSBMLToString(document)
end
