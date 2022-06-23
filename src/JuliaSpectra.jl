module JuliaSpectra

export BasisState, LinearBasisState, SymTopBasisState

# Set up BasisState types
abstract type BasisState end

abstract type LinearBasisState <: BasisState end   

abstract type SymTopBasisState <: BasisState end


# Now include specific molecule models
include("MoleculeParameters.jl")
include("TensorOperators.jl")
include("LinearCaseA.jl")
include("LinearCaseA_Field.jl")
include("LinearCaseB_Field.jl")
include("Linear_BendingMode_Field.jl")
include("SymTopCaseB_Field.jl")
include("EnergyLevels.jl")
include("StateUtils.jl")
include("LineLists.jl")
include("PrintUtils.jl")




end