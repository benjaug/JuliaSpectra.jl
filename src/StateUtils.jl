using Parameters
using LinearAlgebra
using CompositeStructs
using NamedTupleTools
using Graphs: connected_components, SimpleGraph


# Make a Hamiltonian type
@with_kw mutable struct Hamiltonian{T}
    basis::Vector{<:BasisState}
    H_operator::T
    M::Array{Float64, 2} = zeros(Float64, length(basis), length(basis))
    #M::Hermitian{ComplexF64,Array{ComplexF64,2}} = Hermitian(fill(complex(0.0,0.0), length(basis), length(basis)))
end
export Hamiltonian


# Make a "HamiltonianBlock" type
@with_kw mutable struct HamiltonianBlock{T<:BasisState}
    basis::Vector{T}
    H::Array{Float64,2}
    #H::Hermitian{ComplexF64, Array{ComplexF64,2}}
    evals::Vector{Float64}
    evecs::Array{Float64}
end
export HamiltonianBlock



# Matirx element for the energy origin of a state.
function Origin(state::BasisState, state′::BasisState)
    if unpack(state) == unpack(state′)
        ME = 1.0
    else
        ME = 0.0
    end

    return ME
end
export Origin

### FUNCTIONS RELATED TO VIBRONIC MANIFOLDS/BLOCKING BY CONNECTED COMPONENTS
Base.@kwdef struct VibronicManifold{T<:BasisState}
    basis::Vector{T}
    prefactors::Vector{Float64}
    Hterms::Vector{Array{Float64,2}}
end
export VibronicManifold

# function makevibronicmanifold(basis::Vector{<:BasisState}, Hpairs)
#     prefactors = zeros(Float64, length(Hpairs))
#     Hterms = [zeros(Float64,length(basis),length(basis)) for _ = 1:length(Hpairs)]
#     for (i,h) in enumerate(Hpairs)
#         prefactors[i] = h[1]
#         Hterms[i] = build(Hamiltonian(basis=basis, H_operator=h[2]))
#     end

#     VibronicManifold(basis=basis, prefactors=prefactors, Hterms=Hterms)
# end

# These functions are meant to be used to generalize the vibronic manifolds to take in matrices rather than matrix element functions.
function makevibronicmanifold(basis::Vector{<:BasisState}, Hpairs, fpairs=[])
	NH = length(Hpairs)
    Nf = length(fpairs)
    prefactors = ones(Float64, length(Hpairs)) #  the H prefactors will be overwritten and the f prefactor should be 1.
    Hterms = [zeros(Float64,length(basis),length(basis)) for _ = 1:(NH+Nf)]
    for (i,h) in enumerate(Hpairs)
        prefactors[i] = h[1]
        Hterms[i] = build(Hamiltonian(basis=basis, H_operator=h[2]))
    end
	
	for (i,fpair) in enumerate(fpairs)
		for (j,state) in enumerate(basis)
			for (k,state′) in enumerate(basis)
				Hterms[i+NH][j,k] = fpair[2](state,state′,fpair[1]) 
			end
		end
	end

    return VibronicManifold(basis=basis, prefactors=prefactors, Hterms=Hterms)
end
export makevibronicmanifold

# update!(vibronicmanifold, fpairs)
# 	basis = vibronicmanifold.basis
# 	for (i,fpair) in enumerate(fpairs)
# 		for (j,state) in enumerate(basis)
# 			for (k,state') in enumerate(basis)
# 				vibronicmanifold.Hterms[i][j,k] = fpair[2](state,state',fpair[1]) 
# 			end
# 		end
# 	end
# end

function makeTDMvibronicmanifolds(ground::VibronicManifold, excited::VibronicManifold)
    basis_g = ground.basis
    basis_e = excited.basis

    TDMDict = Dict{ Int64, Array{Float64,2} }()
    
    TDM_basis = zeros(Float64, length(basis_g), length(basis_e))
    for p = -1:1
        for (i, state) in enumerate(basis_g)
            for (j, state′) in enumerate(basis_e)
                TDM_basis[i,j] = TDM_E1(state, state′,p)
            end
        end
        TDMDict[p] = TDM_basis
        TDM_basis = zeros(Float64, length(basis_g), length(basis_e))
    end

    return TDMDict
end
export makeTDMvibronicmanifolds

function makeblockedvibronicmanifold(state::VibronicManifold; block_fully = true)
    # if block_fully = false, just put together a single Hamiltonian block.
    # if block_fully = true, block to the maximum extent using connected components.
    if block_fully == false
        basis = state.basis
        Hblock = sum(state.prefactors[j] .* state.Hterms[j] for j in 1:length(state.prefactors))
        evals, evecs = eigen(Hblock)
        return HamiltonianBlock(basis, Hblock, evals, evecs), nothing
    elseif block_fully == true
        state_dict = Dict{Int, HamiltonianBlock}() # initialize Dict to hold output.
        Htemplate = (sum(state.Hterms) .!=0) # "template" Hamiltonian that shows all possible couplings by ignoring prefactors.
        g = SimpleGraph(Htemplate)
        ccs = connected_components(g)
        for (i,block) in enumerate(ccs)
            basis = state.basis[block]
            Hblock = sum(state.prefactors[j] .* state.Hterms[j] for j in 1:length(state.prefactors))[block,block]
            evals, evecs = eigen(Hblock)
            #= could add function to better label the blocks here...  =#
            state_dict[i] = HamiltonianBlock(basis, Hblock, evals, evecs)
        end
        return state_dict, ccs
    end
end
export makeblockedvibronicmanifold

function makeblockedTDMDict(ground_ccs, excited_ccs, TDMvibroDict)

    blockTDMDict = Dict{Tuple{Int,Int,Int},Array{Float64,2}}()

    for p = -1:1
        for (ig,g) in enumerate(ground_ccs)
            for (ie,e) in enumerate(excited_ccs)
                blockTDMDict[(ie,ig,p)] = TDMvibroDict[p][g,e]
            end
        end
    end

    return blockTDMDict
end
export makeblockedTDMDict

### END FUNCTIONS RELATED TO VIBRONIC MANIFOLDS/BLOCKING BY CONNECTED COMPONENTS

function update_matrix!(H::Hamiltonian)
    for (i, state) ∈ enumerate(H.basis)
        for (j, state′) ∈ enumerate(H.basis)
            H.M[i,j] = H.H_operator(state, state′)
        end
    end
    return nothing
end
export update_matrix!

# Function to build the Hamiltonian
function build(H::Hamiltonian; kwargs...)
    basis = H.basis
    H_basis = zeros(Float64, length(basis), length(basis))
    #H_basis = Hermitian(fill(complex(0.0,0.0), length(basis), length(basis)))

    args = [kwarg.second for kwarg in kwargs]
    if !isempty(args)
        H_operator = H.H_operator(args...)
    else
        H_operator = H.H_operator
    end
    for (i, state) in enumerate(basis)
        for (j, state′) in enumerate(basis)
            if i >= j
                H_basis[i,j] = H_operator(state, state′)
            end
#             println(H_basis[i,j])
        end
    end
    H_basis = H_basis + H_basis' - Diagonal(diag(H_basis))
    H.M = H_basis #states * H_basis * states'
#     eigvals, eigvecs = eigen(H.M)
end
export build

# Function to build and diagonalize the Hamiltonian
function solve(H::Hamiltonian; kwargs...)
    basis = H.basis
    H_basis = zeros(Float64, length(basis), length(basis))
    #H_basis = Hermitian(fill(complex(0.0,0.0), length(basis), length(basis)))

    args = [kwarg.second for kwarg in kwargs]
    if !isempty(args)
        H_operator = H.H_operator(args...)
    else
        H_operator = H.H_operator
    end
    for (i, state) in enumerate(basis)
        for (j, state′) in enumerate(basis)
            H_basis[i,j] = H_operator(state, state′)
#             println(H_basis[i,j])
        end
    end
    H.M = H_basis #states * H_basis * states'
#     eigvals, eigvecs = eigen(H.M)
    return eigen(H.M)
end
export solve

# Function to make a dictionary of Hamiltonian blocks.
function makeBlockedState(BasisType, H_operator, QN_bounds, BlockBy, BlockRange)

    state_dict = Dict{Rational, HamiltonianBlock}()

    for blockQN in BlockRange
        QNs_block = merge(QN_bounds, [BlockBy => blockQN])
        basis = enumerate_states(BasisType, QNs_block)

        H = Hamiltonian(basis=basis, H_operator = H_operator)
        evals, evecs = solve(H)

        state_dict[blockQN] = HamiltonianBlock(basis, H.M, evals, evecs)
    end

    return state_dict
end
export makeBlockedState

# Function to make the TDM matrix.
function makeTDMmatrix(excited, ground, p=nothing)
    basis_g = ground.basis
    basis_e = excited.basis

    TDM_basis = zeros(Float64, length(basis_g), length(basis_e))

    for (i, state) in enumerate(basis_g)
        for (j, state′) in enumerate(basis_e)
            TDM_basis[i,j] = TDM_E1(state, state′,p)
        end
    end

    return TDM_basis
end
export makeTDMmatrix

function makeTDMDict(ground, excited)

    TDMDict = Dict{ Tuple{Rational,Rational,Int}, Array{Float64,2} }()
     
    for QN_g in keys(ground)
       for QN_e in keys(excited)
           for p = -1:1
               TDMDict[(QN_g, QN_e, p)] = makeTDMmatrix(excited[QN_g], ground[QN_e],p)
           end
       end
    end

    return TDMDict
end
export makeTDMDict

# Function to make a dictionary of Hamiltonian blocks.
function blockbyconnected(H::Hamiltonian)

    state_dict = Dict{Int, HamiltonianBlock}()

    g = SimpleGraph(H.M)
    ccs = connected_components(g)

    for (i,block) in enumerate(ccs)
        basis = H.basis[block]

        Hblock = Hamiltonian(basis=basis, H_operator = H.H_operator)
        evals, evecs = solve(Hblock)

        # if hasfield(typeof(basis[1]),:F)
        #     state_dict[basis[1].F] = HamiltonianBlock(basis, Hblock.M, evals, evecs)
        # elseif hasfield(typeof(basis[1]),:M)
        #     state_dict[basis[1].M] = HamiltonianBlock(basis, Hblock.M, evals, evecs)
        # else
        #     error("Couldn't get the block label!")
        # end
        state_dict[i] = HamiltonianBlock(basis,Hblock.M, evals,evecs)
    end

    return state_dict
end
export blockbyconnected


# Make a Kronecker Delta function.
δ(x,y) = ==(x,y)
export δ

# Define functions that build up a basis set from constraints on quantum numbers.
function enumerate_states(state_type, QN_bounds)
    states = state_type[]
    QNs = fieldnames(state_type)[1:end-1]
    
    # Define a state with all QN = 0 to get the constraints for the QNs
    η = NamedTuple([QN => 0 for QN in QNs])
    
    enumerate_states(η, states, state_type, QNs, QN_bounds, 1)

    return states
end

function enumerate_states(η, states, state_type, QNs, QN_bounds, idx, max_states=1000)

    if length(states) > max_states
        return false
    end

    if idx == length(QNs) + 1
        new_state = state_type(; η...)
        push!(states, new_state)
        return nothing
    end
    
    # Check if quantum number has been given bounds; else apply constraints
    iterated_QN = QNs[idx]

    QN_constraints = state_type(; η...).constraints
    if iterated_QN ∈ keys(QN_constraints)
        QN_constraint = QN_constraints[iterated_QN]
        QN_constraint_bounds = eval(QN_constraint)
    end
    
    if iterated_QN ∈ keys(QN_bounds)
        bounds = QN_bounds[iterated_QN]
        for i ∈ bounds
            if iterated_QN ∉ keys(QN_constraints) || (i ∈ QN_constraint_bounds)
                η′ = (; η..., iterated_QN => i)
                enumerate_states(η′, states, state_type, QNs, QN_bounds, idx + 1)
            end
        end
    elseif iterated_QN ∈ keys(QN_constraints)
        for i ∈ QN_constraint_bounds
            η′ = (; η..., iterated_QN => i)
            enumerate_states(η′, states, state_type, QNs, QN_bounds, idx + 1)
        end
    else
        enumerate_states(η, states, state_type, QNs, QN_bounds, idx + 1)
    end
    
    return nothing
end
export enumerate_states

# Function to convert Case b to Case a basis.
function convertbasis(state::AbstractCaseB, state′::AbstractCaseA)
    # Hirota 2.3.3
    S,  I,  Λ,  N,  J,  F,  M    = state.S, state.I, state.Λ, state.N, state.J, state.F, state.M
    S′, I′, Λ′, Σ, Ω, J′, F′, M′ = state′.S, state′.I, state′.Λ, state′.Σ, state′.Ω, state′.J, state′.F, state′.M
    return (-1)^(-J+Ω+2S) * sqrt(2N + 1) * wigner3j_(J, N, S, Ω, -Λ, -Σ) * δ(S,S′) * δ(I,I′) * δ(M,M′) * δ(F,F′) * δ(J,J′)
end
convertbasis(state::AbstractCaseA, state′::AbstractCaseB) = convertbasis(state′, state)

# function convertbais(state::UncoupledCaseB, state′::AbstractCaseB)
#     S,  I, Λ, N, J, MN, MS, MI   = state.S, state.I, state.Λ, state.N, state.J, state.MN, state.MS, state.MI
#     S′,  I′,  Λ′,  N′,  J′,  F′,  M′    = state.S, state.I, state.Λ, state.N, state.J, state.F, state.M

#     out= 0.0
#     if M′==(MN+MS+MI) && Λ==Λ′ && N==N′ && S==S′ && I==I′
#         MJ = MN + MS
#         MF = MN+MS+MI
#         out = (-1)^(I-J+MF+S-N+MJ) * sqrt((2F′+1)*(2J′+1)) * wigner3j_(J′,I,F′,MJ,MI,-MF) * wigner3j_(N′,S,J′,MN,MS,-MJ)
#     end
#     return out
# end
# convertbais(state::AbstractCaseA, state′::AbstractUncoupled) = convertbasis(state′, state)
export convertbasis


# Extend the base "+" function such that we can simply add matrix elements before calling them
import Base.+, Base.*, Base.^,Base.-
+(f::Function, g::Function) =
    (args...; kwargs...) -> f(args...; kwargs...) + g(args...; kwargs...)
-(f::Function, g::Function) =
    (args...; kwargs...) -> f(args...; kwargs...) - g(args...; kwargs...)
*(f::Function, g::Function) =
    (args...; kwargs...) -> f(args...; kwargs...) * g(args...; kwargs...)
*(c::Number, f::Function) =
    (args...; kwargs...) -> c * f(args...; kwargs...)
^(f::Function, c::Real) = 
    (args...; kwargs...) -> f(args...; kwargs...)^c

# try adding a method to multiply a tuple of constants times a function that outputs a tuple
*(c::Tuple, f::Function) =
    (args...; kwargs...) -> sum(c[i] * f(args...; kwargs...)[i] for i in 1:length(c))
export +, *, ^,-
