
#= NOTE createReservedWords2.jl also has this function
F: get rid of the trailing newlines, comments, blank lines in a file;
I: file path, read line-by-line;
O: array of strings, each element maps to a line in the file;
B: load the whole file at the same time, may require large memory for large file;
A: read line-by-line, ingnore blank lines and pure comments, throw away in-line comments;
=#
function getRidOfNewlineHashcommentBlankline(inputFilePath::AbstractString)
  f = open(inputFilePath, "r")
  lines = readlines(f)
  newlines = Array{Any, 1}()  # container for results
  for (id, val) in enumerate(lines)
    val = chomp(val)  # trailing newline
    hashId = findfirst(isequal('#'), val)
    if something(hashId,-1) != -1  # get rid of comments
      val = val[1:hashId-1]
    end
    if length(split(val)) != 0  # non-empty
      # throw away ugly trailing spaces, although can be done when splitting
      while(val[end] == ' ')
        val = val[1:end-1]
      end
      push!(newlines, (val, id))
      println(val)
    end
  end
  # println(">"^20)
  # println(newlines)
  return newlines
end


#= NOTE tokenization -> assume can be done just by spaces
FIOBA: using "split" directly
=#
function sentenceTokenization1(inputSenArray::Array)
  tokenizedInputSentencesArray = Array{Array{String,1}, 1}()
  for sen in inputSentencesArray
    tokens = split(sen)
    push!(tokenizedInputSentencesArray, tokens)
  end
  println(typeof(tokenizedInputSentencesArray), size(tokenizedInputSentencesArray))
  return tokenizedInputSentencesArray
end
#= NOTE
F: tokenizing a sentence by processing it char-by-char.
I: output of getRidOfNewlineHashcommentBlankline.
O: Array of Array of tokens.
B: .
A: tokenizing a sentence by processing it char-by-char, punctuation and left/right
   parenthesis are treated as independent tokens, other tokens are seperated by
   space or punctuation set.
=#
function sentenceTokenization2(inputSenArray::Array)
  puncSet = Set(['|', '&', ',', '(', ')'])
  tokenizedInputSentencesArray = Array{Any, 1}()
  for (sen, id) in inputSenArray
    sen = sen*" "  # add a space to end the sentence
    tokens = Array{String,1}()
    currentToken = ""
    for char in sen
      if char == ' '  # encounter " ", push! if have someting
        if !isempty(currentToken)
          push!(tokens, currentToken)
          currentToken = ""
        end
      elseif in(char, puncSet)
        if !isempty(currentToken)
          push!(tokens, currentToken)
          currentToken = ""
        end
        push!(tokens, string(char))
      else
        currentToken *= string(char)
        # println(currentToken)
      end
    end
    push!(tokenizedInputSentencesArray, (tokens, id))
    println(tokens)
  end
  println(typeof(tokenizedInputSentencesArray), size(tokenizedInputSentencesArray))
  return tokenizedInputSentencesArray
end


#= Token classifier 1:
F: get rid of negligible reserved words, tagging or normalizing
I: array{String, 1}, usually from getRidOfNewlineHashcommentBlankline
O: Array of Arrays of Tuples --> Tuple(token, tag)
B: tag all unknown tokens as BioSym, @FIXME how to format for grammar searching?
A: use Dict[tag] = Set{words}, e.g. createReservedWords1.jl
=#
function tokenClassification1(sentencesArray::Array, reservedWords::Dict)
  taggedSentencesArray = Array{Array, 1}()  # container for tagged array
  TAGs = ["induce", "repress", "BioSym"]  # can be incorporated into createReservedWords.jl
  for tokens in sentencesArray  # tag each sentence
    taggedTokens = Array{Tuple,1}()  # Tuple(token, tag)
    for token in tokens  # tag each token in a sentence
      if !(token in reservedWords["negligible"])  # significant words
        breakEndsForLoop = false
        for tag in TAGs  # find a tag
          if token in reservedWords[tag]
            push!(taggedTokens, (token, tag))
            breakEndsForLoop = true
            break
          end
        end
        if !breakEndsForLoop  # token of unknown tagging type FIXME: add a boolean value to indicate status
          # push!(taggedTokens, (token, "UNK"))
          push!(taggedTokens, (token, "BioSym"))  # TODO tag general biological symbols
        end
      end
    end
    push!(taggedSentencesArray, taggedTokens)
  end
  return taggedSentencesArray
end
#= Token classifier 2:
F: get rid of negligible reserved words, tagging or normalizing
I: array{String, 1}, usually from getRidOfNewlineHashcommentBlankline
O: Array of Arrays of Tuples --> Tuple(token, tag)
B: tag all unknown tokens as BioSym, @FIXME how to format for grammar searching?
A: use Dict[word] = tag, e.g. createReservedWords2.jl
=#
function tokenClassification2(sentencesArray::Array, reservedWordsPath::AbstractString)
  include(reservedWordsPath)
  taggedSentencesArray = Array{Any, 1}()  # container for tagged array
  for (tokens, id) in sentencesArray  # tag each sentence
    taggedTokens = Array{Tuple,1}()  # Tuple(token, tag)
    for token in tokens  # tag each token in a sentence
      if haskey(reservedWords, token)  # a reserved word
        tag = reservedWords[token]
        if tag != "negligible"  # important reserved word --> keep it
          push!(taggedTokens, (token, tag))
        end
      else  # non-reserved word --> BioSym
        push!(taggedTokens, (token, "BioSym"))
      end
    end
    push!(taggedSentencesArray, (taggedTokens, id))
  end
  return taggedSentencesArray
end


# # NOTE test
# println("\n---------loading sentences------------")
# inputSentencesArray = getRidOfNewlineHashcommentBlankline("../test/fbacase.txt")
# tokenizedInputSentencesArray = sentenceTokenization2(inputSentencesArray)
#
# println("\n------------normalizing and tagging------------")
# # include("reservedWords.jl")
# taggedSentencesArray = tokenClassification2(tokenizedInputSentencesArray,
#                                             "reservedWords.jl")
# print(taggedSentencesArray)
# println("-----------")
# # reshape for output observation
# printTagSen = [["$(y[1])/$(y[2])" for y in x[1]] for x in taggedSentencesArray]
# println(typeof(printTagSen))
# for ts in printTagSen
#   println(join(ts, "  "))
# end
