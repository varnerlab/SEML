module SEML

include("preprocessor3.jl")
include("biosymDecoder.jl")
include("semanticChecking.jl")
include("modelGeneration.jl")
include("JuliaStrategy.jl")
include("Python2Strategy.jl")
include("Python3Strategy.jl")
include("MATLABStrategy.jl")
include("make_model.jl")

# Export the interface function -
export make_model

end # module
