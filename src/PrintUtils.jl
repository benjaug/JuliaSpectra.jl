using LaTeXStrings

# function print_basis_state(state, fields)
#     str = "\\left|"
#     for (i, field) in enumerate(fields)
#         val = getfield(state, field)
#         val_str = string(Float64(val))
#         str *= string(field) * " = " * val_str 
#         if i < length(fields)
#             str *= string(", \\hspace{0.5em}")
#         end
#     end
#     str *= "\\right\\rangle"
#     return str
# end

# function Base.show(io::IO, m::MIME"text/latex", state::BasisState)

#     basis_type = typeof(state)
#     fields = fieldnames(basis_type)[1:end-1]

#     str = print_basis_state(state, fields)
#     latex_str = latexstring(str)

#     println(io, latex_str)

#     return nothing
# end


function print_basis_state(state, fields)
    str = "|"
    for (i, field) in enumerate(fields)
        val = getfield(state, field)
        val_str = string(Float64(val))
        str *= string(field) * " = " * val_str 
        if i < length(fields)
            str *= string(", ")
        end
    end
    str *= ">"
    return str
end

function Base.show(io::IO,state::BasisState)

    basis_type = typeof(state)
    fields = fieldnames(basis_type)[1:end-1]

    str = print_basis_state(state, fields)

    print(io, str)

    return nothing
end

function Base.show(io::IO, state::Eigenstate)

    energy = state.E
    basis = state.basis
    coeffs = state.coeffs
    basis_type = typeof(basis[1])

    printed = false
    fields = fieldnames(basis_type)[1:end-1] # Christian has [2:end-1]... why?
    print(io,"E = "*string(round(energy,digits=6))*"\n")
    for (i,coeff) in enumerate(coeffs)
        str = ""
        basis_state = basis[i]
        if norm(coeff) > 1e-3
            plus_sign = true
            state_str = ""
            real_val = string(round(real(coeff), digits=4))
            imag_val = string(round(imag(coeff), digits=4))
            if (abs(real(coeff)) > 1e-3) && (abs(imag(coeff)) < 1e-5)
                state_str *= real_val
                if real(coeff) < 0
                    plus_sign = false
                end
            elseif (abs(imag(coeff)) > 1e-3) && abs((real(coeff)) < 1e-5)
                state_str *= imag_val * "i"
                if imag(coeff) < 0
                    plus_sign = false
                end
            else
                state_str *= "(" * real_val * " + " * imag_val * "i" * ")"
            end
            if plus_sign && printed
                str *= "+ "
            end
            str *= state_str
            str *= " "
            str *= print_basis_state(basis_state, fields)
            println(io, str)
            printed = true
        end
    end
    return nothing
end
