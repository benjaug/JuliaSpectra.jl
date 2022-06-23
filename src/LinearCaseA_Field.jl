using Parameters
using CompositeStructs
using InfiniteArrays


export 
    LinearCaseA_Field,
    unpack,
    Rotation, SpinRotation, SpinOrbit, ΛDoubling_p2q, HyperfineIL, HyperfineIS, Hyperfine_IF, Hyperfine_Dipolar_c,
    Zeeman_L, Zeeman_S, Zeeman_parity, Zeeman_gl,
    parity, TDM_E1

Base.@kwdef struct LinearCaseA_Field <: AbstractCaseA #LinearBasisState
    Λ::Rational{Int64}
    S::Rational{Int64}
    Σ::Rational{Int64}
    J::Rational{Int64}
    Ω::Rational{Int64}
    I::Rational{Int64}
    F::Rational{Int64}
    M::Rational{Int64}

    constraints = (
        Σ = -S:S,
        Ω = max(-J, Λ+Σ):min(J, Λ+Σ),
        F = abs(J-I):abs(J+I),
        M = -F:F
    )
end


function unpack(state::LinearCaseA_Field)
    return (state.Λ, state.S, state.Σ, state.J, state.Ω, state.I, state.F, state.M)
end


function Rotation(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    Rotational hamiltonian, includes the diagonal (J^2-Jz^2+S^2-Sz^2) 
    and the spin-uncoupling (-2(JxSx+JySy)) parts.
    Spin-uncoupling part taken from Brown and Carrington Eq. 8.364
    """

    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    Diag = ( J*(J+1) - Ω^2 + S*(S+1) - Σ^2 ) * δ(Λ,Λ′) * δ(Ω,Ω′) * δ(Σ,Σ′) * δ(J,J′) * δ(I,I′) * δ(F,F′) * δ(M,M′)

    if δ(Λ, Λ′) && δ(J, J′) && δ(F, F′) && δ(M, M′)
        SpinUnc = (-1)^(J - Ω + S - Σ) * sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
            sum(
                wigner3j_(J, 1, J′, -Ω, q, Ω′) * 
                wigner3j_(S, 1, S, -Σ, q, Σ′)
                for q in [-1,1]
            ) * δ(Λ, Λ′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
    else
        SpinUnc = 0.0
    end

    ME = Diag - 2*SpinUnc
    return ME 
end

function SpinRotation(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    Spin-rotation Hamiltonian. Includes a diagonal part JzSz-S^2 
    and an off-diagonal part JxSx + JySy.
    """

    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    Diag = (Σ^2 - S*(S+1) ) * δ(Σ,Σ′) * δ(Ω,Ω′) * δ(J,J′) * δ(M,M′) * δ(F,F′) * δ(Λ,Λ′)

    if δ(Λ, Λ′)&& δ(J, J′) && δ(F, F′) && δ(M, M′)
        OffDiag = (-1)^(J - Ω + S - Σ) * sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
        sum(
            wigner3j_(J, 1, J′, -Ω, q, Ω′) * 
            wigner3j_(S, 1, S, -Σ, q, Σ′)
            for q in [-1,1]
        ) * δ(Λ, Λ′) * δ(J, J′) * δ(F, F′) * δ(M, M′)
    else
        OffDiag = 0.0
    end

    ME = Diag + OffDiag
    return ME
end




function SpinOrbit(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    Spin-orbit Hamiltonian, only the diagonal part.
    """

    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    ME = ( Λ*Σ )
    return ME * δ(Λ,Λ′) * δ(Ω,Ω′) * δ(Σ,Σ′) * δ(J,J′) * δ(I,I′) * δ(F,F′) * δ(M,M′)
end


function ΛDoubling_p2q(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    Λ-doubling term parameterized by p+2q. Taken from Brown and Carrington Eq. 9.66
    """

    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(F,F′) && δ(M, M′) && δ(J,J′)
        ME = sqrt( J*(J+1)*(2J+1) * S*(S+1)*(2S+1) ) * 
            sum( 
                δ(Λ′,Λ+2q) * (-1)^(J-Ω+S-Σ) * wigner3j_(J,1,J,-Ω,-q,Ω′) * wigner3j_(S,1,S,-Σ,q,Σ′) 
                for q in [-1,1]
            )
    else
        ME = 0.0
    end

    return ME * δ(F,F′) * δ(M, M′) * δ(J,J′)
end


function HyperfineIL(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    Hyperfine I.L term. Taken from Brown and Carringotn Eq. X.XX
    """

    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(Σ,Σ′) && δ(M,M′) && δ(F,F′) && δ(Ω,Ω′) && δ(Λ,Λ′)
        ME = Λ * (-1)^(J′+I+F+J-Ω) * sqrt(I*(I+1)*(2I+1)*(2J+1)*(2J′+1)) * 
            wigner6j_(J′,I,F,I,J,1) * wigner3j_(J,1,J′,-Ω,0,Ω′)
    else
        ME = 0.0
    end

    return ME * δ(Σ,Σ′) * δ(M,M′) * δ(F,F′) * δ(Ω,Ω′) * δ(Λ,Λ′)
end



function HyperfineIS(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    Hyperfine I.S term. Taken from ...
    """

    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(Λ,Λ′) && δ(F,F′) && δ(M,M′)
        ME = (-1)^(I+J′+F+S-Σ+J-Ω) * sqrt(I*(I+1)*(2I+1)*(2J+1)*(2J′+1)*S*(S+1)*(2S+1)) * 
            wigner6j_(J′,I,F,I,J,1) * 
            sum(
                wigner3j_(J,1,J′,-Ω,q,Ω′) * wigner3j_(S,1,S,-Σ,q,Σ′) 
                for q in -1:1
            )
    else
        ME = 0.0
    end

    return ME * δ(Λ,Λ′) * δ(F,F′) * δ(M,M′)
end


function Hyperfine_IF(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    # Fermi contact interaction
    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(F, F′) && δ(M, M′)
        ME =  (-1)^(I + J′ + F + S - Σ + J - Ω) * 
            sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1)) *
            wigner6j_(I, J′, F, J, I, 1) *
            sum(
                wigner3j_(J, 1, J′, -Ω, q, Ω′) *
                wigner3j_(S, 1, S, -Σ, q, Σ′)
                for q in -1:1
            ) *
            δ(F, F′) * δ(M, M′)
    else
        ME = 0.0
    end

    return ME
end

    
function Hyperfine_Dipolar_c(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(F, F′) && δ(M, M′)
        ME = sqrt(30) * (1/3) * (-1)^(I + J′ + F) * (-1)^(J - Ω) * (-1)^(S - Σ) *
            wigner6j_(I, J′, F, J, I, 1) *
            sqrt( I * (I + 1) * (2I + 1) ) *
            sqrt( (2J + 1) * (2J′ + 1) ) *
            sqrt( S * (S + 1) * (2S + 1) ) *
            sum(
                (-1)^q * 
                wigner3j_(J, 1, J′, -Ω, q, Ω) *
                sum(
                    wigner3j_(1, 2, 1, q′, 0, -q) *
                    wigner3j_(S, 1, S, -Σ, q′, Σ′)
                    for q′ in -1:1
                ) for q in -1:1
            ) * 
            δ(F, F′) * δ(M, M′)
    else
        ME = 0.0
    end

    return ME
end


function Zeeman_L(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    Electron orbit Zeeman.
    """
    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(Ω,Ω′) && δ(Σ,Σ′) && δ(Λ,Λ′) && δ(M,M′) && δ(S,S′) && δ(I,I′)
        ME = Λ * (-1)^(F-M+F′+J+I+1+J-Ω) * wigner6j_(J,F,I,F′,J′,1) * 
            wigner3j_(F,1,F′,-M,0,M) * sqrt((2F+1)*(2F′+1)*(2J+1)*(2J′+1)) * 
            wigner3j_(J,1,J′,-Ω,0,Ω)
    else
        ME = 0.0
    end

    return ME * δ(Ω,Ω′) * δ(Σ,Σ′) * δ(Λ,Λ′) * δ(M,M′) * δ(S,S′) * δ(I,I′)
end


function Zeeman_S(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    Electron spin Zeeman.
    """
    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(S,S′) && δ(I,I′) && δ(M,M′) && δ(Λ,Λ′) 
        ME = sum( (-1)^(F-M+J+I+F′+1+J-Ω+S-Σ) * wigner6j_(J,F,I,F′,J′,1) * 
                wigner3j_(F,1,F′,-M,0,M) * sqrt((2F+1)*(2F′+1)*(2J+1)*(2J′+1)*S*(S+1)*(2S+1)) * 
                wigner3j_(J,1,J′,-Ω,q,Ω′) * wigner3j_(S,1,S,-Σ,q,Σ′)
                for q in -1:1
            )
    else
        ME = 0.0
    end

    return ME * δ(S,S′) * δ(I,I′) * δ(M,M′) * δ(Λ,Λ′) 

end

function Zeeman_parity(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    parity dependent Zeeman. 
    Derived myself, but compare to Brown and Carrington Eq . 9.71 . 
        The extra (-1) prefactor comes from the phase convention for ⟨Λ|e^{∓2iϕ}|Λ⟩.
    """
    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(M,M′) && δ(S,S′) && δ(I,I′)
        ME = (-1)*(-1)^(F-M+S-Σ+J-Ω) * (-1)^(F′+J+I+1) * sqrt(S*(S+1)*(2S+1)) * 
        sqrt((2F+1)*(2F′+1)*(2J+1)*(2J′+1)) * wigner3j_(F,1,F′,-M,0,M′) * 
        wigner6j_(J′,F′,I,F,J,1) * 
            sum(
                δ(Λ′,Λ-2q) * 1 * wigner3j_(S,1,S,-Σ,-q,Σ′) * wigner3j_(J,1,J′,-Ω,q,Ω′)
                for q in [-1,1]
            )
    else
        ME = 0.0
    end

    return ME * δ(M,M′) * δ(S,S′) * δ(I,I′)

end

function Zeeman_gl(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    """
    From CaH optical Zeeman paper, Eq. 4
    """
    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if δ(Λ,Λ′) && δ(S,S′) && δ(I,I′)
        ME = (-1)^(F-M) * wigner3j_(F,1,F′,-M,0,M′) * (-1)^(F′+J+I+1) * sqrt((2F+1)*(2F′+1)) * 
            wigner6j_(J′,F′,I,F,J,1) * 
            sum(
                (-1)^(J-Ω) * wigner3j_(J,1,J′,-Ω,q,Ω′) * sqrt((2J+1)*(2J′+1)) *
                (-1)^(S-Σ) * wigner3j_(S,1,S,-Σ,q,Σ′) * sqrt((2S+1)*(2S′+1))
                for q in [-1,1]
            )
    else
        ME = 0.0
    end

    return ME * δ(Λ,Λ′) * δ(S,S′) * δ(I,I′)

end


function parity(state::LinearCaseA_Field, state′::LinearCaseA_Field)
    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    ME = 0.0
    if Λ′ == -Λ && J'==J && Σ′ == -Σ && Ω′ == -Ω && F==F′ && M==M′
        ME = (-1)^(J-1/2)
    end

    return ME

end


function TDM_E1(state::LinearCaseA_Field, state′::LinearCaseA_Field, p::Int64)
    """
    Transition dipole moment operator (basically Stark operator).
    Polarization has spherical component labeled by p.
    """

    Λ, S, Σ, J, Ω, I, F, M = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′, M′ = unpack(state′)

    if -M′+p+M == 0
        ME = (-1)^p * (-1)^(F′-M′) * wigner3j_(F′,1,F, -M′, p, M) * 
            (-1)^(F+J′+I+1) * sqrt((2F′+1)*(2F+1)) * wigner6j_(J,F,I,F′,J′,1) *
            sum(
                δ(Σ,Σ′) * (-1)^(J′-Ω′) * sqrt((2J′+1)*(2J+1)) * wigner3j_(J′,1,J,-Ω′,q,Ω)
                for q in -1:1 # used to be for q in [-1,1], but this didn't work for Sigma-Sigma transitions.
            )
    else
        ME = 0.0
    end

    return ME
end





