######## Part 1 Useful Information Retrieval ###################

#=
F: extracting BioSym & Verb of each sentence;
I: each sentence in Array of tuples, sentenceType dictionary;
O: sentence type, & array of biosymbol groups;
BA
=#
function extractBioSymGroups(tokenTagPairs::Array, keyVerbs::Set)
  BioSymGroupsArray = []
  typeMarker = ""
  parameterSetting = []
  i = 1
  while i <= length(tokenTagPairs)
    if tokenTagPairs[i][2] == "BioSym"  # find a BioSym group
      tmpBioSymGroup = [tokenTagPairs[i][1]]
      # repeat until find the whole group
      while (i+1 <= length(tokenTagPairs)) && (tokenTagPairs[i+1][2] == "BioSym")
        i += 1  # go for next one
        push!(tmpBioSymGroup, tokenTagPairs[i][1])
      end
      push!(BioSymGroupsArray, tmpBioSymGroup)
      # i += 1  # go for next one
    elseif in(tokenTagPairs[i][2], keyVerbs)
      typeMarker = tokenTagPairs[i][2]
      # i += 1  # go for next one
    elseif tokenTagPairs[i][2] == "reversible"
      push!(parameterSetting, "reversible")
      # i += 1  # go for next one
    elseif tokenTagPairs[i][2] == "irreversible"
      push!(parameterSetting, "irreversible")
    end
    i += 1  # go for next one
  end
  if typeMarker == ""
    println("ERROR: SENTENCE TYPE UNDEFINED for: $(tokenTagPairs)")
  end
  return typeMarker, BioSymGroupsArray, parameterSetting
end

#=
F: extracting BioSym & Key Verb information from tagged sentences, call
   extractBioSymGroups() to process each sentence;
I: Array of sentence, with each sentence in form of array of tuples (token, tag);
O: Array of Tuples, with each tuple in form of (sentence type, array of biosymbol groups);
BA
=#
function extractBioSymVerbInformation(taggedSenteces::Array,
                                      sentenceTypeKeyVerbs::Set, error::Dict)
  BioSymGroupsKeyVerbArray = []
  for i in 1:length(taggedSenteces)
    senType, BioSym, paraSet = extractBioSymGroups(taggedSenteces[i][1], sentenceTypeKeyVerbs)
    lineNum = taggedSenteces[i][2]
    if length(senType) == 0 # catch SENTENCE TYPE UNDEFINED error
        if haskey(error, lineNum)
            push!(error[lineNum], "ERROR: SENTENCE TYPE UNDEFINED")
        else
            error[lineNum] = ["ERROR: SENTENCE TYPE UNDEFINED"]
        end
    else  # successfully extract useful information
        if length(paraSet) == 0
          push!(BioSymGroupsKeyVerbArray, ((senType, lineNum), BioSym))
        else
          push!(BioSymGroupsKeyVerbArray, ((senType, lineNum), BioSym, paraSet))
        end
    end
  end
  return BioSymGroupsKeyVerbArray
end


########## Part 2 Symbol Decoding ###############################

#=
F: decoding biosymbol groups by enumerating all possible combinations
I: Array of tuples (Verb, biosymbol groups), "biosymbol groups" is array of
   biosymbol arrays.
O: biosymbol groups are replaced with corresponding array of biosymbol array of
   biosymbol arrays. (2-layer array for simple case, 3-layer array for symbol
   groups with parenthesis)
BA
=#
function decodingBioSymGroups(VerbBioSymArray::Array, error::Dict)
  # println(VerbBioSymArray)
  isID = Set{Int64}()
  for (id, vbsa) in enumerate(VerbBioSymArray)
    lineNum = vbsa[1][2]
    for i in 1:length(vbsa[2])
      errSym = join(vbsa[2][i], " ")
      vbsa[2][i], err = decodingABioSymGroupByLogicalRules(vbsa[2][i])
      if length(err) != 0
          push!(isID, id)
          if !haskey(error, lineNum)
              error[lineNum] = []
          end
          for e in err
              push!(error[lineNum], "in $(errSym), " * e)
          end
      end
    end
  end
  # delete this from array if has errors?
  deleteat!(VerbBioSymArray, isID)
end

#= ambiguity problem
F: just handling simple cases (one-layer of parenthesis),
   e.g. (A or B) & C & D & (E or F), (A & B) or C or D.
   --> currently, try to do the right thing, but does not consider all possible cases.
I: array of biological tokens (including logcial tokens & bio-tokens)
O: Array of Strings for simple cases; or array of arrays of strings for cases
   with parenthesis.
BA
=#
function decodingABioSymGroupByLogicalRules(tokens::Array)
  println("\n\ninside decodingABioSymGroupByLogicalRules")
  println(tokens)
  err_mes = []  # collection of error messages
  andSet = Set(["and", "&", ","])
  orSet = Set(["or", "|"])

  finalBioSymGroups = []  # elements in parallel relation
  BioSymGroup = []  # for insiders, each element is in the form of [logicalRelation, token1, token2, ...]

  tokensOutParenthesis = []  # for the outsider
  i = 1
  while i <= length(tokens)
    if tokens[i] == "("  # find a ()
      tokensInsideParenthesis = []  # for tokens inside a pair of ()
      i += 1
      while (i<=length(tokens)) && (tokens[i] != ")") # assume one layer of ()
        push!(tokensInsideParenthesis, tokens[i])
        i += 1
      end
      # out of while loop, i is at ")" if i <= length(tokens)
      if i <= length(tokens) && tokens[i] == ")" && length(tokensInsideParenthesis) != 0
        #  has tokens insider () -> handle it
        logicInParenthesis, tmpBioSymGroup, errMes = decodingASequenceOfTokensWithOneLogicalRelation(
        tokensInsideParenthesis, andSet, orSet)
        if errMes != ""  # error occurs
            push!(err_mes, errMes * "inside ()")
        else
          if logicInParenthesis == ""
              logicInParenthesis = "and"
          end
          prepend!(tmpBioSymGroup, [logicInParenthesis])
          push!(BioSymGroup, tmpBioSymGroup)
        end
      else  # what kind of errors?
          push!(err_mes, "ERRORS: incomplete or empty group inside \'()\'")
          println("ERRORS: incomplete or empty group inside \'()\'")
      end
    else  # the outsiders
      push!(tokensOutParenthesis, tokens[i])
    end
    i += 1  # move forward
  end
  # handling outsiders
  logicOutParen, outsiderGroup, errMes = decodingASequenceOfTokensWithOneLogicalRelation(
  tokensOutParenthesis, andSet, orSet)
  if errMes != ""  # error occurs
      push!(err_mes, errMes * "outside ()")
  else
    println("outsiderGroup: $outsiderGroup")
    println("logicOutParen: $logicOutParen")
    if logicOutParen == ""  # no relation for outsiders
      if length(outsiderGroup) == 0  # nothing outside ()
        if length(BioSymGroup) == 1 # can only have 1 group
          finalBioSymGroups = reshapeBioSymbolGroups([], BioSymGroup[1][2:end], BioSymGroup[1][1])
          # println("Simple Bio-sym group")
        else
          println("ERROR: unspecified/unclear relationship")
          push!(err_mes, "ERROR: unspecified/unclear relationship")
        end
      elseif length(outsiderGroup) == 1 && length(BioSymGroup) == 0 # only one token
        finalBioSymGroups = deepcopy(outsiderGroup)
        println("ONE token biosymbol")
      else
        println("ERROR: unspecified/unclear relationship")
        push!(err_mes, "ERROR: unspecified/unclear relationship")
      end
    else  # check logics in/outside () and reshape BioSymGroup for return
      goOnAllow =  true
      for (id, biosym) in enumerate(BioSymGroup)  # check logics in/outside ()
        if biosym[1] == logicOutParen  # logic check
          push!(err_mes, "ERROR: logic inside and outside \'()\' should be different")
          println("ERROR: logic inside and outside \'()\' should be different")
          goOnAllow =  false
        end
      end
      if goOnAllow  # reshape
        BioSymGroup = [subg[2:end] for subg in BioSymGroup]  # reshape
        finalBioSymGroups = reshapeBioSymbolGroups(BioSymGroup, outsiderGroup, logicOutParen)
        println(finalBioSymGroups)
      end
    end
    # println("\nfinalBioSymGroups: ")
    # println(s for s in finalBioSymGroups)
  end

  return finalBioSymGroups, err_mes
end

#=
F: take in a sequence of tokens, seperate bio tokens and logical tokens, return
   logical relation and bio symbol array. Assume only one relation in this sequence.
IOBA
=#
function decodingASequenceOfTokensWithOneLogicalRelation(tokens::Array, andSet::Set, orSet::Set)
  # println("inside decodingASequenceOfTokensWithOneLogicalRelation")
  biosymGroup = []
  logicalRelation = [""]
  errorMessage = ""
  for i = 1:length(tokens)
    if in(tokens[i], andSet)  # find "and"
      if logicalRelation[1] == ""
        logicalRelation[1] = "and"
      elseif logicalRelation[1] == "or"
        errorMessage = "ERROR: LOGICAL INCONSISTENCY"
        println("LOGICAL INCONSISTENCY")
        return logicalRelation[1], BioSymGroup, errorMessage
      else
        continue
      end
    elseif in(tokens[i], orSet)  # find "or"
      if logicalRelation[1] == ""
        logicalRelation[1] = "or"
      elseif logicalRelation[1] == "and"
          errorMessage = "ERROR: LOGICAL INCONSISTENCY"
          println("LOGICAL INCONSISTENCY")
          return logicalRelation[1], BioSymGroup, errorMessage
      else
        continue
      end
    else
      push!(biosymGroup, tokens[i])
    end
  end
  # println(logicalRelation[1])
  # println("jump out decodingASequenceOfTokensWithOneLogicalRelation")
  return logicalRelation[1], biosymGroup, errorMessage
end

#=
F: get all permutations from "OR"/"AND" relations
IOBA
=#
function reshapeBioSymbolGroups(inParenGroups, outParenGroup, logicOutParen)
  finalGroups = []
  if logicOutParen == "and"
    println("inParenGroups: $inParenGroups")
    println("outParenGroup: $outParenGroup")
    if length(inParenGroups) == 0  # simple case: Array of strings
      finalGroups = deepcopy(outParenGroup)
      println("finalGroups: $finalGroups")
      println("Logic AND simple bio-symbol group")
    else  # complicated case: array of arrays of strings
      finalGroups = [outParenGroup]
      # println("finalGroups: $finalGroups")
      for i = 1:length(inParenGroups)
        tmp2 = []
        for biosym in inParenGroups[i]  # append each token to the previous finalGroups
          tmp1 = deepcopy(finalGroups)
          # println("tmp1: $tmp1")
          for j = 1:length(tmp1)
            push!(tmp1[j], biosym)
          end
          # println("tmp1: $tmp1")
          append!(tmp2, tmp1)
        end
        # println("tmp2: $tmp2")
        finalGroups = deepcopy(tmp2)  # update finalGroups after processing each outParenGroup element
        # println("finalGroups: $finalGroups")
      end
      println("Logic AND complicated bio-symbol group")
    end
  else  # "or" --> split outsiders into arrays, add insiders
    for biosym in outParenGroup
      push!(finalGroups, [biosym])
    end
    if length(inParenGroups) == 0 # simple case
      # finalGroups = outParenGroup
      println("Logic OR simple bio-symbol group")
    else
      for biosymgroup in inParenGroups
        push!(finalGroups, biosymgroup)
      end
      println("Logic OR complicated bio-symbol group")
    end
  end
  return finalGroups
end

#=
F: print out an array of tuples, each tuple is a pair of a token and a multiple-
   layer arrays.
I: an array of tuples, each tuple is a pair of a token and a multiple-
   layer arrays. Array of arrays of strings (2-layers) for simple cases,
   array of arrays of arrays of strings (3-layers) for complicated cases.
O: see terminal.
BA
=#
function printArrayOfTripleToken(tokenPlusMultiLayerArray::Array, targetDataType::Type)
  println("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  for tuple in tokenPlusMultiLayerArray
    println(tuple[1])
    println()
    if length(tuple[2]) == 0  # empty biosymbols
      println("THIS sentence is incorrect")
    else  # non-empty
      for symArray1 in tuple[2]  # 1st layer
        if length(symArray1) == 0
          println("THIS sentence is incorrect")
        else  # FIXME: handle default return data structure '[]'
          if typeof(symArray1[1]) <: targetDataType || length(symArray1[1]) == 0
            # println("print 1: simple case --> Array of strings")
            println(symArray1)
          else
            for symArray2 in symArray1  # 2nd layer
              if typeof(symArray2[1]) <: targetDataType
                # println("print 2: complicated case --> Array of arrays")
                println(symArray2)
              else
                for symArray3 in symArray2  # 3rd layer
                  if typeof(symArray3[1]) <: targetDataType
                    println("print 3: too ambiguous to handle")
                    println(symArray3)
                  else
                    for symArray4 in symArray3  # 4th layer
                      if typeof(symArray4[1]) <: targetDataType
                        println("print 4: too ambiguous to handle")
                        println(symArray4)
                      else
                        println(symArray4)
                        println("NEED PRINT 5: world of no man")
                      end
                    end
                  end
                end
              end
            end
          end
        end
        println("")
      end
    end
    if isdefined(tuple, 3)
      println(tuple[3])
    end
    println("-----------------")
  end
end


########### Part 3 Symbol Conversion Dictionary & general biosym replacement ####

#=
F: set up symbol conversion table, delete symbol conversion sentences
I: verb-SymbolGroup pairs
O: conversion dictionary, also modify the input array;
BA
=#
function setUpSymbolConvertionDict(verbSymArray::Array, error::Dict)
  # verbSymArray = Array of pairs of ((senType, line#), BioSym, paraSet)
  conversionTable = Dict{AbstractString, AbstractString}()
  isID = Array{Int64, 1}()
  for (id, vs) in enumerate(verbSymArray)
    if vs[1][1] == "is"  #  find conversion sentence --> check legality
      if ((length(vs[2]) == 2) && (typeof(vs[2][1]) <: Array) && (typeof(vs[2][2]) <: Array)
        && (typeof(vs[2][1][1]) <: AbstractString) && (typeof(vs[2][2][1]) <: AbstractString)
        && (length(vs[2][1]) == length(vs[2][2])))
        for i = 1:length(vs[2][1])
          conversionTable[vs[2][1][i]] = vs[2][2][i]
        end
      else  # incorrect
        lineNum = vs[1][2]
        if !haskey(error, lineNum)
          error[lineNum] = []
        end
        push!(error[lineNum], "incorrect symbol conversion sentence")
        print("incorrect symbol conversion sentence found")
      end
      push!(isID, id)
    end
  end
  deleteat!(verbSymArray, isID)
  conversionTable["[]"] = "[]"
  conversionTable["SITE"] = "SITE"
  return conversionTable
end

mutable struct generalBioSym
  oriBioName::AbstractString
  bioType::AbstractString
  coeff::Float64  # Float64 has default value 0.0?
  function generalBioSym()
    new()
  end
end

#=
F: decompose a token into generalBioSym, check coefficient and type
I: a token, conversion table
O: a generalBioSym if success; otherwise, an empty array []
BA
=#
function BioSym2GeneralBioSym(token::AbstractString, convDict::Dict, err::Array)
  bioSym = generalBioSym()
  if token == "[]"  # special symbol '[]'
    bioSym.oriBioName = "[]"
    bioSym.bioType = "[]"
    bioSym.coeff = 0.0
  elseif (!occursin("*", token) && !occursin("_", token))  # site
    bioSym.oriBioName = token
    bioSym.bioType = "SITE"
    bioSym.coeff = 1.0
  else  # normal symbol
    tokenPiece1 = split(token, "*")
    if length(tokenPiece1) == 0 || length(tokenPiece1) > 2  # "*" is reserved for specifying coefficient
      println("error while parsing coefficient for \'$token\'")
      push!(err, "error while parsing coefficient for \'$token\'")
    else  # maybe 1) correct 2) default 3) error
      if length(tokenPiece1) == 2
        if isnothing(tryparse(Float64, tokenPiece1[1]))  # parse() failed
          println("error while parsing coefficient for \'$token\'")
          push!(err, "error while parsing coefficient for \'$token\'")
        else  # got name and coefficient
          bioSym.coeff = parse(Float64, tokenPiece1[1])
          bioSym.oriBioName = tokenPiece1[2]
        end
      else  # default coefficient
        bioSym.coeff = 1.0
        bioSym.oriBioName = tokenPiece1[1]
      end
      # how to continue?
      if isdefined(bioSym, :oriBioName)  # means success in finding coefficient
        tokenPiece2 = split(bioSym.oriBioName, "_")
        if length(tokenPiece2) < 3  # no type or no name specified
          println("no type/compartment specified for \'$token\'")
          push!(err, "no type/compartment specified for \'$token\'")
        else
          tmpType = tokenPiece2[1]*"_"
          if tmpType in keys(convDict)  # want to change into system-defined type? -> YES
            bioSym.bioType = tmpType
          else
            println("incorrect type for \'$token\'")
            push!(err, "incorrect type for \'$token\'")
          end
        end
      end
    end
  end
  # only return when all fields are parsed correctly, otherwise, return empty array []
  if (isdefined(bioSym, :oriBioName) && isdefined(bioSym, :bioType)) # && bioSym.coeff != 0
    return bioSym
  else
    return []
  end
end

#=
F: replace bio-symbol in input array with general biosymbol
I: array of biosymbols or array of arrays of biosymbols, conversion table
O: array of generlBiosymbols or array of arrays of generalBiosymbols
BA
=#
function BioSymArray2generalBioSymArray(inputArr::Array, convDict::Dict, errM::Array)
  postInputArr = []
  if typeof(inputArr[1]) <: AbstractString  # one-layer (simple case)
    postInputArr = [BioSym2GeneralBioSym(s, convDict, errM) for s in inputArr]
  else
    for arr in inputArr
      if typeof(arr[1]) <: AbstractString  # two-layer (complicated case)
        push!(postInputArr, [BioSym2GeneralBioSym(s, convDict, errM) for s in arr])
      else  # three-layer or more, no consider
        println("Currently in land of no man")
        push!(errM, "sorry, currently in land of no man")
      end
    end
  end
  return postInputArr
end

#=
F: replace bio-symbol in input array with general biosymbol
I: array of verb-biosymbolArray pairs, conversion table
O: array of verb-generalBioSymArray pairs
BA
=#
function replaceEachBioSymWithGeneralBioSym(biosymArray::Array,
                     conversionDict::Dict, error_report::Dict)
  # biosymArray = Array of pairs of ((senType, line#), BioSym, paraSet)
  newVerbSymArray = []
  for bs in biosymArray
    err_mes = []
    lineNum = bs[1][2]
    if ((bs[1][1] == "is") || (bs[1][1] == ""))  # fail to delete it
      println("failed to modify the input array?")
      push!(err_mes, "failed to modify the input array in previous steps")
    else
      tmpArr = [BioSymArray2generalBioSymArray(s, conversionDict, err_mes) for s in bs[2]]
      if length(err_mes) == 0
        if isdefined(bs, 3)  # has parameter setting
          push!(newVerbSymArray, (bs[1], tmpArr, bs[3]))
        else
          push!(newVerbSymArray, (bs[1], tmpArr))
        end
      end
    end
    # error occurs
    if length(err_mes) != 0
      if !haskey(error_report, lineNum)
        error_report[lineNum] = []
      end
      append!(error_report[lineNum], err_mes)
    end
  end
  return newVerbSymArray
end



# # NOTE test
# include("preprocessor3.jl")
# println("\n---------loading sentences------------")
# inputSentencesArray = getRidOfNewlineHashcommentBlankline("../test/fbacase.txt")
# tokenizedInputSentencesArray = sentenceTokenization2(inputSentencesArray)
# println("\n------------normalizing and tagging------------")
# reservedWordsPath = "reservedWords.jl"
# taggedSentencesArray = tokenClassification2(tokenizedInputSentencesArray, reservedWordsPath)
# # reshape for output observation
# printTagSen = [["$(y[1])/$(y[2])" for y in x[1]] for x in taggedSentencesArray]
# println(typeof(printTagSen))
# foreach(println, [join(ts, "  ") for ts in printTagSen])
#
# println("\n------------information extraction------------")
# include(reservedWordsPath)
# error_report = Dict{Int64, Array}()
# BioSymVerbInfo = extractBioSymVerbInformation(taggedSentencesArray, reservedWords["SentenceType"], error_report)
# # println("print by foreach:")
# foreach(println, BioSymVerbInfo)
# println("\n------------decoding bio symbols------------")
# decodingBioSymGroups(BioSymVerbInfo, error_report)
# println("\n------------decoding bio symbols results------------")
# printArrayOfTripleToken(BioSymVerbInfo, AbstractString)
# println("\n---------type conversion dictionary---------------")
# typeConversionDict = setUpSymbolConvertionDict(BioSymVerbInfo, error_report)
# foreach(println, typeConversionDict)
# println("\n------replace each biosymbol string with generalBioSym composite-----------------")
# newVerbBiosymInfo = replaceEachBioSymWithGeneralBioSym(BioSymVerbInfo, typeConversionDict, error_report)
# printArrayOfTripleToken(newVerbBiosymInfo, generalBioSym)
