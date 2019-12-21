#= semantic checking
e.g. mRNA_F activates g_H → mRNA can not activate the expression of genes
→ be careful! Biology is error-prone, never be too sure.
=#

#=
F: 1) check biological meaning of each token; 2) have all components for each type,
  get ready for IR generation;
I: array of triplets (verb - generalBioSym Arrays - parameters),
   conversion table, sentence specification (encoded in code);
O: can get rid of type and parameters? why shall?
   Array of tuples, tuple = (reaction type, [reactants, products, catalysts]).
   [reactants, products, catalysts] are all in List form
BA
=#
function semanticCheckingForEachTriplet(tripletArr::Array, conversionDict::Dict,
         error::Dict)
  # tripletArr = Array of ((senType, line#), BioSym, paraSet)
  preIRForm = []  # array of tuple ((verb, line#), reactant, product, catalyst)
  sys2userDict = Dict()
  for (key, val) in conversionDict
    sys2userDict[val] = key
  end

  for tp in tripletArr
    err_mess = []  # collection of error messages
    lineNum = tp[1][2]
    tmpIR = []
    if tp[1][1] == "react"
      tmpIR = checkReactType(tp, conversionDict, err_mess)
    elseif tp[1][1] == "uptake"
      tmpIR = checkUptakeType(tp, conversionDict, sys2userDict, err_mess)
    elseif tp[1][1] == "secrete"
      tmpIR = checkSecreteType(tp, conversionDict, sys2userDict, err_mess)
    elseif tp[1][1] == "bind"
      tmpIR = checkBindType(tp, conversionDict, err_mess)
    elseif tp[1][1] == "unbind"
      tmpIR = checkUnbindType(tp, conversionDict, err_mess)
    elseif tp[1][1] == "phosphorylate"
      tmpIR = checkPhosphorylateType(tp, conversionDict, sys2userDict, err_mess)
    elseif tp[1][1] == "dephosphorylate"
      tmpIR = checkDephosphorylateType(tp, conversionDict, sys2userDict, err_mess)
    elseif tp[1][1] == "catalyze"
      tmpIR = checkCatalyzeType(tp, conversionDict, err_mess)
    elseif tp[1][1] == "induce"  # GENE
      tmpIR = checkInduceType(tp, conversionDict, err_mess)
    elseif tp[1][1] == "activate"  # PROTEIN
      tmpIR = checkInduceType(tp, conversionDict, err_mess)
    elseif tp[1][1] == "repress"  # GENE
      tmpIR = checkRepressType(tp, conversionDict, err_mess)
    elseif tp[1][1] == "inhibit"  # PROTEIN
      tmpIR = checkRepressType(tp, conversionDict, err_mess)
    else
      println("**********unknown sentence type: $(tp[1])***************")
      push!(err_mess, "**********unknown sentence type: $(tp[1])***************")
    end
    if length(err_mess) != 0  # error occur
      if !haskey(error, lineNum)
        error[lineNum] = []
      end
      append!(error[lineNum], err_mess)
    else
      if length(tmpIR) != 0  # have something to return
        append!(preIRForm, tmpIR)
      end
    end
  end
  return preIRForm
end


#==========Functions For Each Type========================#

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkDephosphorylateType(senInfo::Tuple, u2sDict::Dict, s2uDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "[]"])
  legalTypeSet2 = Set(["PROTEIN"])
  legalTypeSet3 = Set(["SITE"])
  # check parameter setting
  rever = false
  if isdefined(senInfo, 3)  # has parameters, may have model parameters in future
    if senInfo[3][1] == "reversible"
      rever =  true
    end
  end
  if DS10SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2) # correct reactant
    if isassigned(senInfo[2], 3)  # have site
      if DS10SemanticChecking(senInfo[2][3], u2sDict, legalTypeSet3)  # correct SITE
        product = createProductsList(senInfo[2][2], "-"*senInfo[2][3][1].oriBioName, "", err)
      else  # incorrect site
        println("*semantic error or incorrect symbol format found in $(senInfo[2][3])")
        push!(err, "*semantic error or incorrect symbol format found in $(senInfo[2][3])")
        return []
      end
    else  # does not have site --> default product
      product = createProductsList(senInfo[2][2], s2uDict["-pho"], "", err)
    end
    if DS1SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1) # correct catalyst of DS1
      if rever  # reversible
        return [(senInfo[1], [senInfo[2][2], product, senInfo[2][1]]), (senInfo[1], [product, senInfo[2][2], senInfo[2][1]])]
      else  # irreversible
        return [(senInfo[1], [senInfo[2][2], product, senInfo[2][1]])]
      end
    elseif DS2SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1) # correct catalyst of DS2
      if rever  # reversible
        returnArr = []
        for cata in senInfo[2][1]
          append!(returnArr, [(senInfo[1], [senInfo[2][2], product, senInfo[2][1]]),
            (senInfo[1], [product, senInfo[2][2], senInfo[2][1]])])
        end
        return returnArr
      else  # irreversible
        returnArr = [(senInfo[1], [senInfo[2][2], product, cata]) for cata in senInfo[2][1]]
        return returnArr
      end
    else  # incorrect catalyst
      println("**semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
      push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
    end
  else  # incorrect reactant
    println("***semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
    push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkPhosphorylateType(senInfo::Tuple, u2sDict::Dict, s2uDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "[]"])
  legalTypeSet2 = Set(["PROTEIN"])
  legalTypeSet3 = Set(["SITE"])
  # check parameter setting
  rever = false
  if isdefined(senInfo, 3)  # has parameters, may have model parameters in future
    if senInfo[3][1] == "reversible"
      rever =  true
    end
  end
  if DS10SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2) # correct reactant
    if isassigned(senInfo[2], 3)  # have site
      if DS10SemanticChecking(senInfo[2][3], u2sDict, legalTypeSet3)  # correct SITE
        product = createProductsList(senInfo[2][2], "-"*senInfo[2][3][1].oriBioName)
      else  # incorrect site
        println("*semantic error or incorrect symbol format found in $(senInfo[2][3])")
        push!(err, "*semantic error or incorrect symbol format found in $(senInfo[2][3])")
        return []
      end
    else  # does not have site --> default product
      product = createProductsList(senInfo[2][2], s2uDict["-pho"])
    end
    if DS1SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1) # correct catalyst of DS1
      if rever  # reversible
        return [(senInfo[1], [senInfo[2][2], product, senInfo[2][1]]), (senInfo[1], [product, senInfo[2][2], senInfo[2][1]])]
      else  # irreversible
        return [(senInfo[1], [senInfo[2][2], product, senInfo[2][1]])]
      end
    elseif DS2SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1) # correct catalyst of DS2
      if rever  # reversible
        returnArr = []
        for cata in senInfo[2][1]
          append!(returnArr, [(senInfo[1], [senInfo[2][2], product, senInfo[2][1]]),
            (senInfo[1], [product, senInfo[2][2], senInfo[2][1]])])
        end
        return returnArr
      else  # irreversible
        returnArr = [(senInfo[1], [senInfo[2][2], product, cata]) for cata in senInfo[2][1]]
        return returnArr
      end
    else  # incorrect catalyst
      println("**semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
      push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
    end
  else  # incorrect reactant
    println("***semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
    push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkCatalyzeType(senInfo::Tuple, u2sDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "[]"])
  legalTypeSet2 = Set(["METABOLITE", "PROTEIN", "GENE"])
  # check parameter setting
  rever = false
  if isdefined(senInfo, 3)  # has parameters, may have model parameters in future
    if senInfo[3][1] == "reversible"
      rever =  true
    end
  end
  if !isassigned(senInfo[2],2) || !isassigned(senInfo[2], 3)
    push!(err, "no enough inputs for \'catalyze\' type")
  else
    if (DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2) &&
      DS1SemanticChecking(senInfo[2][3], u2sDict, legalTypeSet2))  # reactant/product correct
      if DS1SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # catalyst correct
        if rever  # reversible
          returnArr = [(senInfo[1], [senInfo[2][2], senInfo[2][3], senInfo[2][1]]),
            (senInfo[1], [senInfo[2][3], senInfo[2][2], senInfo[2][1]])]
          return returnArr
        else  # irreversible
          returnArr = [(senInfo[1], [senInfo[2][2], senInfo[2][3], senInfo[2][1]])]
          return returnArr
        end
      elseif DS2SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # catalyst correct
        if rever  # reversible
          returnArr = []
          for cata in senInfo[2][1]
            append!(returnArr, [(senInfo[1], [senInfo[2][2], senInfo[2][3], cata]),
              (senInfo[1], [senInfo[2][3], senInfo[2][2], cata])])
          end
          return returnArr
        else  # irreversible
          returnArr = [(senInfo[1], [senInfo[2][2], senInfo[2][3], cata]) for cata in senInfo[2][1]]
          return returnArr
        end
      else  # incorrect catalyst
        println("*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
        push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
      end
    else  # incorrect reactants and/or products
      if !DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)
        println("**semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
        push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
      end
      if !DS1SemanticChecking(senInfo[2][3], u2sDict, legalTypeSet2)
        println("***semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][3])")
        push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][3])")
      end
    end
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkRepressType(senInfo::Tuple, u2sDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "METABOLITE"])
  # Induce or Activate
  if senInfo[1][1] == "repress"
    legalTypeSet2 = Set(["GENE"])
  else
    legalTypeSet2 = Set(["PROTEIN"])
  end
  if DS1SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # regulators correct
    if DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)  # targets correct
      returnArr = [(senInfo[1], [senInfo[2][1], [tar]]) for tar in senInfo[2][2]]
      return returnArr
    else  # semantic error or incorrect targets
      println("*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
      push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
    end
  elseif DS2SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # regulators correct
    if DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)  # targets correct
      returnArr = [(senInfo[1], [regu, [tar]]) for regu in senInfo[2][1] for tar in senInfo[2][2]]
      return returnArr
    else  # semantic error or incorrect targets
      println("**semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
      push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
    end
  else  # semantic error or incorrect targets
    println("***semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
    push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkInduceType(senInfo::Tuple, u2sDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "METABOLITE"])
  # Induce or Activate
  if senInfo[1][1] == "induce"
    legalTypeSet2 = Set(["GENE"])
  else
    legalTypeSet2 = Set(["PROTEIN"])
  end
  if DS1SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # regulators correct
    if DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)  # targets correct
      returnArr = [(senInfo[1], [senInfo[2][1], [tar]]) for tar in senInfo[2][2]]
      return returnArr
    else  # semantic error or incorrect targets
      println("*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
      push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
    end
  elseif DS2SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # regulators correct
    if DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)  # targets correct
      returnArr = [(senInfo[1], [regu, [tar]]) for regu in senInfo[2][1] for tar in senInfo[2][2]]
      return returnArr
    else  # semantic error or incorrect targets
      println("**semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
      push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
    end
  else  # semantic error or incorrect targets
    println("***semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
    push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkUnbindType(senInfo::Tuple, u2sDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "METABOLITE", "MRNA"])
  legalTypeSet2 = Set(["PROTEIN", "METABOLITE", "MRNA", "GENE"])
  # check parameter setting
  rever = true
  if isdefined(senInfo, 3)  # has parameters, may have model parameters in future
    if senInfo[3][1] == "irreversible"
      rever =  false
    end
  end
  if DS10SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # reactants correct
    if isassigned(senInfo[2], 2)  # have products
      if DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)  # correct products
        if rever  # reversible reaction
          returnArr = [(senInfo[1], [senInfo[2][1], senInfo[2][2]]), (senInfo[1], [senInfo[2][2], senInfo[2][1]])]
          return returnArr
        else  # irreversible
          returnArr = [(senInfo[1], [senInfo[2][1], senInfo[2][2]])]
          return returnArr
        end
      else  # incorrect product
        println("*semantic error or depreciated expression in $(senInfo[2][2])")
        push!(err, "*semantic error or depreciated expression in $(senInfo[2][2])")
      end
    else  # error: product required
      println("*for \'unbind\', products need to be specified")
      push!(err, "*for \'unbind\', products need to be specified")
    end
  else  # incorrect reactants
    println("**semantic error or incorrect symbol format found in $(senInfo[2][1])")
    push!(err, "*semantic error or incorrect symbol format found in $(senInfo[2][1])")
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkBindType(senInfo::Tuple, u2sDict::Dict, err::Array)
  # semantic checking, can not all be GENE
  semantic1 = false
  if typeof(senInfo[2][1]) <: Array && typeof(senInfo[2][1][1]) <: generalBioSym
    for gb in senInfo[2][1]
      if u2sDict[gb.bioType] != "GENE"
        semantic1 = true
        break
      end
    end
  end
  legalTypeSet2 = Set(["PROTEIN", "METABOLITE", "MRNA"])
  # check parameter setting
  rever = true
  if isdefined(senInfo, 3)  # has parameters, may have model parameters in future
    if senInfo[3][1] == "irreversible"
      rever =  false
    end
  end
  if semantic1  # reactants correct
    if isassigned(senInfo[2], 2)  # have products
      if DS10SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)  # correct product
        if rever  # reversible reaction
          returnArr = [(senInfo[1], [senInfo[2][1], senInfo[2][2]]), (senInfo[1], [senInfo[2][2], senInfo[2][1]])]
          return returnArr
        else  # irreversible
          returnArr = [(senInfo[1], [senInfo[2][1], senInfo[2][2]])]
          return returnArr
        end
      else  # incorrect product
        println("*semantic error or depreciated expression in $(senInfo[2][2])")
        push!(err, "*semantic error or depreciated expression in $(senInfo[2][2])")
      end
    else  # generate product
      product = createBindProductsList(senInfo[2][1])
      if rever  # reversible reaction
        returnArr = [(senInfo[1], [senInfo[2][1], product]), (senInfo[1], [product, senInfo[2][1]])]
        return returnArr
      else  # irreversible
        returnArr = [(senInfo[1], [senInfo[2][1], product])]
        return returnArr
      end
    end
  else  # incorrect reactants
    println("**semantic error or incorrect symbol format found in $(senInfo[2][1])")
    push!(err, "*semantic error or incorrect symbol format found in $(senInfo[2][1])")
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkSecreteType(senInfo::Tuple, u2sDict::Dict, s2uDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "METABOLITE"])
  legalTypeSet2 = Set(["PROTEIN"])
  # check parameter setting
  rever = false
  if isdefined(senInfo, 3)  # has parameters, may have model parameters in future
    if senInfo[3][1] == "reversible"
      rever =  true
    end
  end
  # semantics = DS1SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)
  if DS1SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # reactants correct
    products = createProductsList(senInfo[2][1], s2uDict["_c"], s2uDict["_exc"], err)
    if isassigned(senInfo[2], 2)  # have transporters
      if DS10SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # unambiguous
        # semantics2 = DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)
        if DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)  # transporters correct
          if rever  # reversible reaction
            returnArr = []
            for s in senInfo[2][2]
              append!(returnArr, [(senInfo[1], [senInfo[2][1], products, [s]]),
                (senInfo[1], [products, senInfo[2][1], [s]])])
            end
            return returnArr
          else  # irreversible
            returnArr = []
            for s in senInfo[2][2]
              append!(returnArr, [(senInfo[1], [senInfo[2][1], products, [s]])])
            end
            return returnArr
          end
        else  # transporter(s) not in the required form
          println("*semantic error or incorrect symbol format found in $(senInfo[2][2])")
          push!(err, "*semantic error or incorrect symbol format found in $(senInfo[2][2])")
        end
      else  # ambiguous
        println("**depreciated expression due to inherent ambiguity in $(senInfo[2][2])")
        push!(err, "*depreciated expression due to inherent ambiguity in $(senInfo[2][2])")
      end
    else  # does not have transporters
      returnArr = []
      if rever  # reversible
        for (id, sym) in enumerate(senInfo[2][1])
          append!(returnArr, [(senInfo[1], [[sym], [products[id]]]), (senInfo[1], [[products[id]], [sym]])])
        end
      else  # irreversible
        for (id, sym) in enumerate(senInfo[2][1])
          append!(returnArr, [(senInfo[1], [[sym], [products[id]]])])
        end
      end
      return returnArr
    end
  else  # incorrect reactants
    println("***semantic error or incorrect symbol format found in $(senInfo[2][1])")
    push!(err, "*semantic error or incorrect symbol format found in $(senInfo[2][1])")
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table, reverse conversion table;
O: array of bi-tuple;
BA
=#
function checkUptakeType(senInfo::Tuple, u2sDict::Dict, s2uDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "METABOLITE", "MRNA", "GENE"])
  legalTypeSet2 = Set(["PROTEIN"])
  # check parameter setting
  rever = false
  if isdefined(senInfo, 3)  # has parameters, may have model parameters in future
    if senInfo[3][1] == "reversible"
      rever =  true
    end
  end
  semantics = DS1SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)
  if semantics  # reactants correct
    products = createProductsList(senInfo[2][1], s2uDict["_exc"], s2uDict["_c"], err)
    if isassigned(senInfo[2], 2)  # have transporters
      if DS10SemanticChecking(senInfo[2][1], u2sDict, legalTypeSet1)  # unambiguous
        semantics2 = DS1SemanticChecking(senInfo[2][2], u2sDict, legalTypeSet2)
        if semantics2  # transporters correct
          if rever  # reversible reaction
            returnArr = []
            for s in senInfo[2][2]
              append!(returnArr, [(senInfo[1], [senInfo[2][1], products, [s]]),
                (senInfo[1], [products, senInfo[2][1], [s]])])
            end
            return returnArr
          else  # irreversible
            returnArr = []
            for s in senInfo[2][2]
              append!(returnArr, [(senInfo[1], [senInfo[2][1], products, [s]])])
            end
            return returnArr
          end
        else  # transporter(s) not in the required form
          println("*semantic error or incorrect symbol format found in $(senInfo[2][2])")
          push!(err, "*semantic error or incorrect symbol format found in $(senInfo[2][2])")
        end
      else  # ambiguous
        println("**depreciated expression due to inherent ambiguity in $(senInfo[2][2])")
        push!(err, "*depreciated expression due to inherent ambiguity in $(senInfo[2][2])")
      end
    else  # does not have transporters
      returnArr = []
      if rever  # reversible
        for (id, sym) in enumerate(senInfo[2][1])
          append!(returnArr, [(senInfo[1], [[sym], [products[id]]]), (senInfo[1], [[products[id]], [sym]])])
        end
      else  # irreversible
        for (id, sym) in enumerate(senInfo[2][1])
          append!(returnArr, [(senInfo[1], [[sym], [products[id]]])])
        end
      end
      return returnArr
    end
  else  # incorrect reactants
    println("***semantic error or incorrect symbol format found in $(senInfo[2][1])")
    push!(err, "*semantic error or incorrect symbol format found in $(senInfo[2][1])")
  end
  return []
end

#=
F:
I: Tuple of biosymbols and parameters, conversion table;
O: array of bi-tuple;
BA
=#
function checkReactType(senInfo::Tuple, conDict::Dict, err::Array)
  legalTypeSet1 = Set(["PROTEIN", "METABOLITE", "MRNA", "GENE"])
  legalTypeSet2 = Set(["[]"])
  # check parameter setting
  rever = true
  if isdefined(senInfo, 3)  # has parameters, may have model parameters in future
    if senInfo[3][1] == "irreversible"
      rever =  false
    end
  end
  # must have symbols
  if length(senInfo[2]) < 2
    println("no enough input for sentence $(senInfo[1][1])")
    push!(err, "no enough input for sentence $(senInfo[1][1])")
  else
    if length(senInfo[2]) > 2
      println("too many inputs, only the first two groups are handled")
      push!(err, "too many inputs, only the first two groups are handled")
    end
    if DS10SemanticChecking(senInfo[2][1], conDict, legalTypeSet2)  # comes from nothing
      if DS1SemanticChecking(senInfo[2][2], conDict, legalTypeSet1)  # products correct
        if rever  # reversible
          returnArr = []
          for reac in senInfo[2][2]
            append!(returnArr, [(senInfo[1], [senInfo[2][1], [reac]]),
              (senInfo[1], [[reac], senInfo[2][1]])])
          end
          return returnArr
        else  # irreversible
          returnArr = [(senInfo[1], [senInfo[2][1], [reac]]) for reac in senInfo[2][2]]
          return returnArr
        end
      else  # product incorrect
        println("*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
        push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][2])")
      end
    elseif DS10SemanticChecking(senInfo[2][2], conDict, legalTypeSet2)  # goes to nothing
      if DS1SemanticChecking(senInfo[2][1], conDict, legalTypeSet1)  # products correct
        if rever  # reversible
          returnArr = []
          for reac in senInfo[2][1]
            append!(returnArr, [(senInfo[1], [[reac], senInfo[2][2]]),
              (senInfo[1], [senInfo[2][2], [reac]])])
          end
          return returnArr
        else  # irreversible
          returnArr = [(senInfo[1], [[reac], senInfo[2][2]]) for reac in senInfo[2][1]]
          return returnArr
        end
      else  # product incorrect
        println("**semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
        push!(err, "*semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1])")
      end
    else  # incorrect usage of "react" or error
      println("***incorrect use of \'react\', suggest using \'catalyze\';" *
        " OR semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1]) and $(senInfo[2][2])")
      push!(err, "***incorrect use of \'react\', suggest using \'catalyze\';" *
        " OR semantic error or ambiguity or incorrect symbol format found in $(senInfo[2][1]) and $(senInfo[2][2])")
    end
  end
  return []
end


#==========Create Products List=============================#

#=
F: take in an array, create binding compounds, link by hyphen
IOBA
=#
function createBindProductsList(oriList::Array)
  product = generalBioSym()
  product.coeff = 1.0
  product.bioType = oriList[1].bioType  # FIXME compound type -> mRNA-pro, DNA-pro
  product.oriBioName = join([s.oriBioName for s in oriList], "--")
  return [product]
end

#=
F: take in an array, create corresponding products list
IOBA
=#
function createProductsList(oriList::Array, oldTag, newTag, err)
  newList = []
  for arr in oriList
    new_arr = deepcopy(arr)
    if !occursin(oldTag, new_arr.oriBioName)
      println("ERROR: no \'$oldTag\' found in \'$(new_arr.oriBioName)\'")
      push!(err, "ERROR: no \'$oldTag\' found in \'$(new_arr.oriBioName)\'")
    end
    new_arr.oriBioName = replace(new_arr.oriBioName, oldTag => newTag)
    push!(newList, new_arr)
  end
  return newList
end

#=
F: take in an array, create corresponding products list
IOBA
=#
function createProductsList(oriList::Array, addTag::AbstractString)
  newList = []
  for arr in oriList
    new_arr = deepcopy(arr)
    new_arr.oriBioName = new_arr.oriBioName * addTag
    push!(newList, new_arr)
  end
  return newList
end

#==========Data Structure Checking=========================#

#=
F: check an input array is an array with a single/multiple element(s) which
  is/are of generalBioSym type
I:
O: semanticStatus: true if correct, otherwise, false;
BA
=#
function DS1SemanticChecking(inputArr::Array, typeDict::Dict, typeSet)
  if (DS10SemanticChecking(inputArr::Array, typeDict::Dict, typeSet) ||
    DS11SemanticChecking(inputArr::Array, typeDict::Dict, typeSet))
    return true
  else
    return false
  end
end

#=
F: check an input array is an array with a single element which is of generalBioSym type
I:
O: semanticStatus: true for correct, false for incorrect;
BA
=#
function DS10SemanticChecking(inputArr::Array, typeDict::Dict, typeSet)
  semanticStatus = true
  if typeof(inputArr) <: Array && length(inputArr) == 1 && typeof(inputArr[1]) <: generalBioSym
    for gb in inputArr
      if !in(typeDict[gb.bioType], typeSet)
        semanticStatus = false
        # println("semantic error found in $gb")
        break
      end
    end
  else  # incorrect format
    semanticStatus = false
    # println("incorrect symbol format for \'react\' type")
  end
  return semanticStatus
end

#=
F: check an input array is an array with many a single element which is of generalBioSym type
I:
O: semanticStatus: true for correct, false for incorrect;
BA
=#
function DS11SemanticChecking(inputArr::Array, typeDict::Dict, typeSet)
  semanticStatus = true
  if typeof(inputArr) <: Array && length(inputArr) > 1 && typeof(inputArr[1]) <: generalBioSym
    for gb in inputArr
      if !in(typeDict[gb.bioType], typeSet)
        semanticStatus = false
        # println("semantic error found in $gb")
        break
      end
    end
  else  # incorrect format
    semanticStatus = false
    # println("incorrect symbol format for \'react\' type")
  end
  return semanticStatus
end

#=
F: check an input array is an array of arrays with a single/multiple element(s)
   which is/are of generalBioSym type.
I:
O: semanticStatus: true if correct, otherwise, false;
BA
=#
function DS2SemanticChecking(inputArr::Array, typeDict::Dict, typeSet)
  semantic = false
  if (typeof(inputArr) <: Array) && (typeof(inputArr[1]) <: Array)  # data structure 2
    semantic = minimum([DS1SemanticChecking(s, typeDict, typeSet) for s in inputArr])
  end
  return semantic
end
