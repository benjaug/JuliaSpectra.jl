using Parameters
using LinearAlgebra
using CompositeStructs
using NamedTupleTools
using LoopVectorization
using PrettyTables
using Plots

export 
    Transition, makeLineList, plotLineList, makeTransitionTable,
    filter_transitions, label_transitions!,
    LineList_tabular

Base.@kwdef mutable struct Transition
    lower::Eigenstate
    upper::Eigenstate
    intensity::Float64
    frequency::Float64
    label::String
    function Transition(lower::Eigenstate, upper::Eigenstate, intensity::Float64, frequency::Float64)
        return new(lower, upper, intensity, frequency, "")
    end
end

# Function to make a LineList from a ground/excited state pair.
function makeLineList(ground, excited, p=nothing; TK=10, min_I = 1e-10)
    LineList = Transition[]
 
    QN_g_min = minimum(keys(ground))
    QN_g_max = maximum(keys(ground))
    QN_e_min = minimum(keys(excited))
    QN_e_max = maximum(keys(excited))

    Elowest = minimum(minimum.([ground[Fg].evals for Fg in QN_g_min:QN_g_max]))
     
    for QN_g = QN_g_min:QN_g_max
        for QN_e = QN_e_min:QN_e_max
            #@time begin
            if p == "par"
                TDM_mat = makeTDMmatrix(excited[QN_e], ground[QN_g], 0)
            elseif p == "perp"
                TDM_mat = 1/sqrt(2)*( makeTDMmatrix(excited[QN_e], ground[QN_g], -1) - makeTDMmatrix(excited[QN_e], ground[QN_g], 1))
            else
                TDM_mat = makeTDMmatrix(excited[QN_e], ground[QN_g],p)
            end
            #end
            size_g_state = length(ground[QN_g].basis)
            size_e_state = length(excited[QN_e].basis)

            for m = 1:size_g_state
                g_vec = ground[QN_g].evecs[:,m]
                g_val = ground[QN_g].evals[m]
                for n = 1:size_e_state
                    e_vec = excited[QN_e].evecs[:,n]
                    e_val = excited[QN_e].evals[n]

                    BF = exp(-(g_val - Elowest)/(TK) * 0.695) # just set temperaature to 10 K and convert to cm-1. Eventually make this a parameter.
                    #line_strength = abs(g_vec' * TDM_mat * e_vec)^2
                    line_strength = abs2(g_vec' * TDM_mat * e_vec)
                    intensity = line_strength * BF

                    if intensity > min_I
                        eigl = Eigenstate(g_val, g_vec, ground[QN_g].basis)
                        eigu = Eigenstate(e_val, e_vec, excited[QN_e].basis)
                        push!(LineList, Transition(eigl, eigu, intensity, e_val-g_val))
                    end
                end
            end
        end
    end

    return LineList
end

# # This version can make a linelist from a precomputed dictionary of TDM elements between basis states.
# function makeLineList(ground, excited, TDMDict::Dict, p=nothing; TK=10, min_I = 1e-10, return_struct = false)
#     if return_struct == false
#         LineList = zeros(Float64, 0, 2)
#     else
#         LineList = Transition[]
#     end
 
#     Elowest = minimum(minimum.([ground[Fg].evals for Fg in keys(ground)]))
     
#     for QN_g in keys(ground)
#         for QN_e = keys(excited)
#             if typeof(p) <: Vector
#                 TDM_mat = sum( p[i] * (-1)^(i-2) * TDMDict[QN_e, QN_g, i-2] for i in 1:3) # order basis as [sigma-, pi, sigma+]. i-2 gives [-1,0,+1] to label these.
#             else    
#                 TDM_mat = TDMDict[QN_e, QN_g, p]
#             end 
#             size_g_state = length(ground[QN_g].basis)
#             size_e_state = length(excited[QN_e].basis)

#             for m = 1:size_g_state
#                 g_vec = ground[QN_g].evecs[:,m]
#                 g_val = ground[QN_g].evals[m]
#                 for n = 1:size_e_state
#                     e_vec = excited[QN_e].evecs[:,n]
#                     e_val = excited[QN_e].evals[n]

#                     BF = exp(-(g_val - Elowest)/(TK) * 0.695) # just set temperaature to 10 K and convert to cm-1. Eventually make this a parameter.

#                     # line_strength = abs(g_vec' * TDM_mat * e_vec)^2
#                     line_strength = abs2(g_vec' * TDM_mat * e_vec)
#                     intensity = line_strength * BF

#                     if intensity > min_I
#                         if return_struct == false
#                             LineList = vcat(LineList, [e_val - g_val intensity])
#                         else
#                             eigl = Eigenstate(g_val, g_vec, ground[QN_g].basis)
#                             eigu = Eigenstate(e_val, e_vec, excited[QN_e].basis)
#                             push!(LineList, Transition(eigl, eigu, intensity, e_val-g_val))
#                         end
#                     end
#                 end
#             end
#         end
#     end

#     return LineList
# end
# This version can make a linelist from a precomputed dictionary of TDM elements between basis states.
function makeLineList(ground, excited, TDMDict::Dict, p=nothing; TK=10, min_I = 1e-10)
    LineList = Transition[]
    Elowest = minimum(minimum.([ground[Fg].evals for Fg in keys(ground)]))
    for QN_g in keys(ground)
        for QN_e = keys(excited)
            if typeof(p) <: Vector
                TDM_mat = sum( p[i] * (-1)^(i-2) * TDMDict[QN_e, QN_g, i-2] for i in 1:3) # order basis as [sigma-, pi, sigma+]. i-2 gives [-1,0,+1] to label these.
            else    
                TDM_mat = TDMDict[QN_e, QN_g, p]
            end 
            size_g_state = length(ground[QN_g].basis)
            size_e_state = length(excited[QN_e].basis)

            for m = 1:size_g_state
                g_vec = ground[QN_g].evecs[:,m]
                g_val = ground[QN_g].evals[m]
                for n = 1:size_e_state
                    e_vec = excited[QN_e].evecs[:,n]
                    e_val = excited[QN_e].evals[n]

                    BF = exp(-(g_val - Elowest)/(TK) * 0.695) # 0.695 converts K to cm-1. 

                    # line_strength = abs(g_vec' * TDM_mat * e_vec)^2
                    line_strength = abs2(g_vec' * TDM_mat * e_vec)
                    intensity = line_strength * BF

                    if intensity > min_I
                        eigl = Eigenstate(g_val, g_vec, ground[QN_g].basis)
                        eigu = Eigenstate(e_val, e_vec, excited[QN_e].basis)
                        push!(LineList, Transition(eigl, eigu, intensity, e_val-g_val))
                    end
                end
            end
        end
    end
    return LineList
end

function filter_transitions(LineList, filter_rule::Function)
    LineList_filtered = Transition[]
    for (_,line) in enumerate(LineList)
        res = filter_rule(line)
        if res == true 
            push!(LineList_filtered, line)
        end
    end
    return LineList_filtered
end

function label_transitions!(LineList, label_rule::Function)
    # LineList = Vector{Transition}
    # label_rule = function that takes a Transition and returns a label
    for (_,trans) in enumerate(LineList)
        trans.label = label_rule(trans)
    end
end

# # Function to plot a LineList 
# function plotLineList(LineList, fmin, fmax; gamma = 0.001, nstep=100_000, interactive_labels = false)
#     step = (fmax-fmin)/(nstep)
#     freq = zeros(Float64,nstep)
#     amp = zeros(Float64,nstep)
#     if interactive_labels == false
#         @turbo for i = 1:size(LineList,1)
#             for j = 1:nstep
#                 offset = j*step
#                 dif = (offset - (LineList[i,1]-fmin))^2
#                 freq[j] = offset + fmin
#                 amp[j] = amp[j] + (gamma/(2*3.141)) * (1/(dif + (gamma^2)/4)) * LineList[i,2]
#             end
#         end
#         return freq, amp
#     elseif interactive_labels == true
#         xd = Float64[]
#         yd = Float64[]
#         text = String[]    #@turbo for (i,line) = enumerate(LineList)
#         for (i,line) = enumerate(LineList)
#             for j = 1:nstep
#                 offset = j*step
#                 dif = (offset - (line.frequency-fmin))^2
#                 freq[j] = offset + fmin
#                 amp[j] = amp[j] + (gamma/(2*3.141)) * (1/(dif + (gamma^2)/4)) * line.intensity
#             end
    
#             push!(xd,line.frequency)
#             push!(yd,line.intensity)
#             push!(text, line.label)
#         end
#         ph = plot(xd, yd, linealpha=0, xlabel="Frequency (MHz)", hover = text, frame=:box)
#         plot!(ph,freq,amp, label="Spectrum", hover = false, legend=false, xlims=(fmin,fmax))
#         return ph
#     end
# end
# Function to plot a LineList 
function plotLineList(LineList, fmin, fmax; gamma = 0.001, nstep=10_000, make_plot = false, kwargs...)
    step = (fmax-fmin)/(nstep)
    freq = zeros(Float64,nstep)
    amp = zeros(Float64,nstep)
    xd = Float64[]
    yd = Float64[]
    text = String[]   
    # Loop through frequency array and make the Lorentzian peaks at each position. 
    # Add all possible transitions together.  
    for (_,line) = enumerate(LineList)
        @turbo for j = 1:nstep
            offset = j*step
            dif = (offset - (line.frequency-fmin))^2
            freq[j] = offset + fmin
            amp[j] = amp[j] + (gamma/(2*3.141)) * (1/(dif + (gamma^2)/4)) * line.intensity
        end
    end
    # now make the list of labels and positions
    for (_,line) = enumerate(LineList)
        push!(xd,line.frequency)
        ind = findnearest(freq, line.frequency)
        push!(yd, amp[ind])
        push!(text, string(round(line.frequency,digits=3),"\n",line.label))
    end
    # different plot output options.
    if make_plot == false
        return freq, amp
    elseif make_plot == true
        ph = Plots.plot(xd, yd, linealpha=0, hover = text, frame=:box)
        Plots.plot!(ph, freq,amp; hover=false, kwargs...)
        return ph
    end
end
# function to find nearest point in array
#useful for labeling the spectrum by the peak that corresponds to the transition underlying it.
findnearest(A::AbstractArray,t) = findmin(abs.(A .- t))[2]

# Function to make a LineList from a ground/excited state pair.
function makeTransitionTable(ground, excited, p=nothing)
 
    QN_g_min = minimum(keys(ground))
    QN_g_max = maximum(keys(ground))
    QN_e_min = minimum(keys(excited))
    QN_e_max = maximum(keys(excited))
     
    for QN_g = QN_g_min:QN_g_max
        for QN_e = QN_e_min:QN_e_max
            TDM_mat = TDMDict[(QN_e, QN_g, p)]
            #@time begin
            if p == "par"
                TDM_mat = makeTDMmatrix(excited[QN_e], ground[QN_g], 0)
            elseif p == "perp"
                TDM_mat = 1/sqrt(2)*( makeTDMmatrix(excited[QN_e], ground[QN_g], -1) - makeTDMmatrix(excited[QN_e], ground[QN_g], 1))
            else
                TDM_mat = makeTDMmatrix(excited[QN_e], ground[QN_g],p)
            end
        end
    end

    return TDM_mat
end

# Output a nicely formatted table generated from the linelist. 
function LineList_tabular(LineList::Vector{Transition}; offset = 0.0, kwargs...)
    sa = StructArray(LineList)
    inds = sortperm(sa.frequency)
    header = (["Transition","Frequency","Intensity"],[" "," - $(offset) (MHz)","(arb.)"])
    pt = pretty_table(hcat(sa[inds].label, (sa[inds].frequency .- offset), sa[inds].intensity); header=header, formatters = ft_round([1,3],[2,3]),
                alignment = :c)
    return pt
end