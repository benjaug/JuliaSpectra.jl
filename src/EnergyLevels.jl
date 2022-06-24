using LinearAlgebra
using Plots
using StructArrays

export 
    eigenshuffle, flatten,
    plottransitiontable, 
    makestateindsdict, 
    transitionintensities, branchingratios,
    subspaceinds, 
    packeigensystem, unpackeigensystem,
    Eigenstate, changebasis,
    plotZeemanlevels, plotZeemanlevels!

struct Eigenstate{T<:BasisState}
    E::Float64
    coeffs::Array{Float64}
    basis::Vector{T}
end
EigenstateStructArray = StructArray{Eigenstate, 2, NamedTuple{(:E, :coeffs, :basis), Tuple{Matrix{Float64}, Matrix{Array{Float64}}, Matrix{Vector{<:BasisState}}}}, Int64}


# Some tools for working with energy levels and making plots.
"""
Sort lists of eigenvectors such that the order preserves maximum wavefunction overlap between the two lists.
"""
function shufflestates(oldvecs, newvecs, r=4; warn_if_lost=true)
    overlap = abs.(round.(oldvecs' * newvecs,digits=r))
    vecdist = sqrt.(abs.(1 .-overlap.^2))
    ordering = [argmin(vecdist[i,:]) for i in 1:size(oldvecs,1)]

    # check we didn't lose any eigenvectors
    if warn_if_lost
        if length(unique(ordering)) !== size(newvecs,1)
            println("Lost an eigenvector!")
        end
    end
    
    return ordering
end

"""
Perform a parameter scan while maintaining consistent eigenvector ordering.
"""
function eigenshuffle(H0, H1, basis, params; pack_eigensystem=true, warn_if_lost=true)
    # scan over a range of parameters while keeping the ordering of eigenvectors consistent.
    # H0 is the "static" part of the Hamiltonian and H1 is the part of the Hamiltonian that is scanned.
    # H1 does NOT include whatever prefactor is being scanned-- that prefactor goes in the list params.
    # The return_raw optional argument decides whether to output the raw eigenvalues/eigenvectors as matrices
    # or whether to return an Eigensystem object.

    basissize = length(basis)
    nsteps = length(params)

    evals = zeros(basissize, nsteps)
    evecs = zeros(basissize, basissize, nsteps)

    Ham0 = Hamiltonian(basis=basis, H_operator = H0)
    build(Ham0)
    Ham1 = Hamiltonian(basis=basis, H_operator = H1)
    build(Ham1)
    H0M = Ham0.M
    H1M = Ham1.M 
    evals[:,1], evecs[:,:,1] = eigen(H0M+params[1]*H1M)

    for i = 2:nsteps
        vals, vecs = eigen(H0M + params[i]*H1M)
        order = shufflestates(evecs[:,:,i-1], vecs, warn_if_lost=warn_if_lost)
        evals[:,i] = vals[order]
        evecs[:,:,i] = vecs[:,order]
    end

    if pack_eigensystem == false
        return evals, evecs
    elseif pack_eigensystem == true
        return packeigensystem(evecs,evals,basis)
    end
end

"""
Take a matrix of eigenvectors and list of eigenvalues and pack it into a list of Eigenstate objects.
"""
function packeigensystem(evecs, evals, basis)
    # Function to build a StructArray holding the Eigenstates.
    nsteps = size(evecs,3)
    nbasis = size(evecs,2)
    s = Array{Eigenstate}(undef,nbasis,nsteps)
    for i = 1:nbasis
        for j = 1:nsteps
            s[i,j] = Eigenstate(evals[i,j], evecs[:,i,j], basis)
        end
    end
    sa = StructArray(s)
    return sa
end

"""
Useful method to turn a vector of eigenvector arrays into a matrix.
"""
flatten(v) = reduce(hcat,v)



function plottransitiontable(strengths, gbasis, ebasis; title=nothing)
    heatmap(strengths', c=palette(:linear_wyor_100_45_c55_n256,rev=false), colorbar=false, legend=:false, framestyle=:box)
    if title !== nothing
        plot!(title=title)
    end

    gstatelines = cumsum((2*1/2+1)*(2*1/2+1)*(2 .* sort(unique(gbasis[i].N for i in 1:length(gbasis))) .+1))
    estatelines = (2*1/2+1)*(2 .* sort(unique(ebasis[i].J for i in 1:length(ebasis))) .+1)
    estatelines = cumsum(repeat(estatelines, outer = [1,2])'[:])
    vline!(0.5 .+ gstatelines, color=:black) # There are (2I+1)(2S+1)(2N+1) sublevels in X.
    hline!(0.5 .+ estatelines, color=:black) # There are (2I+1)(2J+1) sublevels in A, for each parity.
    xlims!(0, length(gbasis)+1)
    ylims!(0, length(ebasis)+1)
    xlabel!("Ground state eigenvalue #")
    ylabel!("Excited state eigenvalue #")
end



function branchingratios(vecsX::Array{Float64}, vecsA::Array{Float64}, TDMDict)
    return abs2.(vecsX' * sum(TDMDict[p] for p in -1:1) * vecsA)
end

function branchingratios(eigsys1::T, eigsys2::T, TDMDict) where T <: StructArray
    if typeof(eigsys1[1].basis) !== typeof(eigsys2[1].basis)
        error("Eigensystem 1 has basis ", typeof(eigsys1[1].basis[1]), " but Eigensystem 2 has basis ", typeof(eigsys2[1].basis[1]))
    end
    vecs1 = flatten(eigsys1.coeffs)
    vecs2 = flatten(eigsys2.coeffs)
    return abs2.(vecs1' * sum(TDMDict[p] for p in -1:1) * vecs2)
end


function transitionintensities(vecsX::Array{Float64}, vecsA::Array{Float64}, TDMDict, p::String)
    if p == "x"
        return abs2.(vecsX' * (TDMDict[1] + TDMDict[-1])* vecsA)
    elseif p == "z"
        return abs2.(vecsX' * (TDMDict[0]) * vecsA)
    end
end

function transitionintensities(vecsX::Array{Float64}, vecsA::Array{Float64}, TDMDict, theta::Number)

    TDMmat = -sin(theta)*(TDMDict[1] + TDMDict[-1]) + cos(theta)*TDMDict[0] # I think a - sign in front of the perpedicular term due to (-1)^p?
    return abs2.(vecsX' * TDMmat * vecsA)

end

function transitionintensities(eigsys1, eigsys2, TDMDict, theta::Float64)
    if typeof(eigsys1[1].basis) !== typeof(eigsys2[1].basis)
        error("Eigensystem 1 has basis ", typeof(eigsys1[1].basis[1]), " but Eigensystem 2 has basis ", typeof(eigsys2[1].basis[1]))
    end
    vecs1 = flatten(eigsys1.coeffs)
    vecs2 = flatten(eigsys2.coeffs)
    TDMmat = -sin(theta)*(TDMDict[1] + TDMDict[-1]) + cos(theta)*TDMDict[0] # I think a - sign in front of the perpedicular term due to (-1)^p?
    return abs2.(vecs1' * TDMmat * vecs2)
end

function subspaceinds(evecs, evals, basis::Vector{<:BasisState}, conditions)
# Return a dictionary that contains indices to the eigenvectors satisfying conditions.
# conditions must have signature: (evecs, evals, basis).
    out = Dict()
    for i = 1:size(evecs,1)
        label = conditions(evecs[:,i,:],evals[i,:],basis)
        if !haskey(out, label)
            out[label] = [i]
        else
            push!(out[label], i)
        end
    end
    return out
end

function subspaceinds(eigsys, conditions) 
    # Return a dictionary that contains indices to the eigenvectors satisfying conditions.
    # conditions must have signature: (StructArray{Eigenstate})
    out = Dict()
    for i = 1:size(eigsys,1)
        label = conditions(eigsys[i,:])
        if !haskey(out, label)
            out[label] = [i]
        else
            push!(out[label], i)
        end
    end
    return out
end

function plotZeemanlevels(eigsys; units="MHz", energy_offset=0.0, kwargs...)
    p = plot(frame=:box, grid=false; kwargs...)
    xlabel!(p,"Magnetic quantum number, M")
    for i in 1:length(eigsys)
        energy = eigsys[i].E
        M = eigsys[i].basis[argmax(abs2.(eigsys[i].coeffs))].M 
        if units == "MHz"
            energy = energy * sol
            ylabel!(p,"Energy (MHz)")
        else 
            ylabel!(p, "Energy (cm-1)")
        end
        plot!(p,[M-0.35, M+0.35],[energy-energy_offset, energy-energy_offset]; kwargs...)
    end
    return p 
end

function plotZeemanlevels!(p,eigsys; units="MHz", energy_offset=0.0, kwargs...)
    p = plot(p,frame=:box, grid=false; kwargs...)
    xlabel!(p,"Magnetic quantum number, M")
    for i in 1:length(eigsys)
        energy = eigsys[i].E
        M = eigsys[i].basis[argmax(abs2.(eigsys[i].coeffs))].M 
        if units == "MHz"
            energy = energy * sol
            ylabel!(p,"Energy (MHz)")
        else 
            ylabel!(p, "Energy (cm-1)")
        end
        plot!(p,[M-0.35, M+0.35],[energy-energy_offset, energy-energy_offset]; kwargs...)
    end
    return p 
end

"""
Change the basis of a StructArray of Eigenstates. 
"""
function changebasis(eigsys, newbasis)
    np = size(eigsys,2) # number of parameter values at which eigensystem was computed.
    oldbasis = eigsys[1,1].basis 

    # Now define the unitary transformation from old to new basis.
    Umat = zeros(length(oldbasis), length(newbasis))
    for (i,state) in enumerate(oldbasis)
        for (j,state′) in enumerate(newbasis)
            Umat[i,j] = convertbasis(state′, state)
        end
    end

    # convert the eigenvectors, looping over parameter values
    evecs = zeros(Float64, length(newbasis), length(oldbasis), np)
    evals = zeros(Float64, length(oldbasis), np)
    for i = 1:np
        evecs[:,:,i] = Umat'*flatten(eigsys[:,i].coeffs) # multiply old basis set by transformation matrix
        evals[:,i] = eigsys[:,i].E # leave energies unchanged
    end

    return packeigensystem(evecs,evals,newbasis)
end