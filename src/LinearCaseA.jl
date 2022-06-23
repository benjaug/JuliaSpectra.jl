using Parameters
using CompositeStructs
using InfiniteArrays

abstract type AbstractCaseA <: LinearBasisState end

Base.@kwdef struct LinearCaseA <: AbstractCaseA
    Λ::Rational{Int64}
    S::Rational{Int64}
    Σ::Rational{Int64}
    J::Rational{Int64}
    Ω::Rational{Int64}
    I::Rational{Int64}
    F::Rational{Int64}

    constraints = (
        Σ = -S:S,
        Ω = max(-J, Λ+Σ):min(J, Λ+Σ),
        # J = abs(F-I):abs(F+I),
        F = abs(J-I):abs(J+I)
    )
end
export LinearCaseA

function unpack(state::LinearCaseA)
    return (state.Λ, state.S, state.Σ, state.J, state.Ω, state.I, state.F)
end
export unpack


function Rotation(state::LinearCaseA, state′::LinearCaseA)
    """
    Rotational hamiltonian, includes the diagonal (J^2-Jz^2+S^2-Sz^2) 
    and the spin-uncoupling (-2(JxSx+JySy)) parts.
    Spin-uncoupling part taken from Brown and Carrington Eq. 8.364
    """

    Λ, S, Σ, J, Ω, I, F = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′ = unpack(state′)

    Diag = 0.0
    if δ(Λ,Λ′) && δ(Ω,Ω′)  && δ(Σ,Σ′) && δ(J,J′) && δ(I,I′) && δ(F,F′)
        Diag = 1.0*( J*(J+1) - Ω^2 + S*(S+1) - Σ^2 )
    end

    SpinUnc = 0.0
    if δ(Λ,Λ′) && δ(J,J′) && δ(F,F′)
        SpinUnc = 1.0*(-1)^(J - Ω + S - Σ) * sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
            sum(
                wigner3j_(J, 1, J′, -Ω, q, Ω′) * 
                wigner3j_(S, 1, S, -Σ, q, Σ′)
                for q in [-1,1]
            ) 
    end

    ME = Diag - 2*SpinUnc
    return ME 
end
export Rotation

function SpinRotation(state::LinearCaseA, state′::LinearCaseA)
    """
    Spin-rotation Hamiltonian. Includes a diagonal part JzSz-S^2 
    and an off-diagonal part JxSx + JySy.
    """

    Λ, S, Σ, J, Ω, I, F = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′ = unpack(state′)

    Diag = 0.0
    if δ(Σ,Σ′) && δ(Ω,Ω′) && δ(J,J′) && δ(F,F′) && δ(Λ,Λ′)
        Diag = (Σ^2 - S*(S+1) )
    end

    OffDiag = 0.0
    if Λ==Λ′ && J==J′ && F==F′
        OffDiag = (-1)^(J - Ω + S - Σ) * sqrt( J * (J + 1) * (2J + 1) * S * (S + 1) * (2S + 1) ) *
                sum(
                    wigner3j_(J, 1, J′, -Ω, q, Ω′) * 
                    wigner3j_(S, 1, S, -Σ, q, Σ′)
                    for q in [-1,1]
                ) * δ(Λ, Λ′) * δ(J, J′) * δ(F, F′)
    end

    ME = Diag + OffDiag
    return ME
end
export SpinRotation



function SpinOrbit(state::LinearCaseA, state′::LinearCaseA)
    """
    Spin-orbit Hamiltonian, only the diagonal part.
    """

    Λ, S, Σ, J, Ω, I, F = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′ = unpack(state′)

    ME = 0.0
    if Λ==Λ′ && Ω==Ω′ && Σ==Σ′ && J==J′ && I==I′ && F==F′
        ME = 1.0*( Λ*Σ )
    end

    return ME 
end
export SpinOrbit

function ΛDoubling_p2q(state::LinearCaseA, state′::LinearCaseA)
    """
    Λ-doubling term parameterized by p+2q. Taken from Brown and Carrington Eq. 9.66
    """

    Λ, S, Σ, J, Ω, I, F = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′ = unpack(state′)


    ME = 0.0
    if F==F′ && J==J′
        ME = sqrt( J*(J+1)*(2J+1) * S*(S+1)*(2S+1) ) * 
            sum( 
                δ(Λ′,Λ+2q) * (-1)^(J-Ω+S-Σ) * wigner3j_(J,1,J,-Ω,-q,Ω′) * wigner3j_(S,1,S,-Σ,q,Σ′) 
                for q in [-1,1]
            )
    end

    return ME
end
export ΛDoubling_p2q

function Hyperfine_IF(state::LinearCaseA, state′::LinearCaseA)
    # Fermi contact interaction
    Λ, S, Σ, J, Ω, I, F= unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′= unpack(state′)

    if δ(F, F′) 
        ME =  (-1)^(I + J′ + F + S - Σ + J - Ω) * 
            sqrt(I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1)) *
            wigner6j_(I, J′, F, J, I, 1) *
            sum(
                wigner3j_(J, 1, J′, -Ω, q, Ω′) *
                wigner3j_(S, 1, S, -Σ, q, Σ′)
                for q in -1:1
            ) *
            δ(F, F′)
    else
        ME = 0.0
    end

    return ME
end
export Hyperfine_IF
    
function Hyperfine_Dipolar_c(state::LinearCaseA, state′::LinearCaseA)
    Λ, S, Σ, J, Ω, I, F= unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′= unpack(state′)

    if δ(F, F′) 
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
            δ(F, F′)
    else
        ME = 0.0
    end

    return ME
end
export Hyperfine_Dipolar_c


function HyperfineIL(state::LinearCaseA, state′::LinearCaseA)
    """
    Hyperfine I.L term. Taken from Brown and Carringotn Eq. X.XX
    """

    Λ, S, Σ, J, Ω, I, F = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′ = unpack(state′)

    ME = 0.0
    if Σ==Σ′ && F==F′ && Ω==Ω′ && Λ==Λ′
        ME = Λ * (-1)^(J′+I+F+J-Ω) * sqrt(I*(I+1)*(2I+1)*(2J+1)*(2J′+1)) * 
            wigner6j_(J′,I,F,I,J,1) * wigner3j_(J,1,J′,-Ω,0,Ω′)
    end

    return ME
end
export HyperfineIL



function HyperfineIS(state::LinearCaseA, state′::LinearCaseA)
    """
    Hyperfine I.S term. Taken from ...
    """

    Λ, S, Σ, J, Ω, I, F = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′ = unpack(state′)

    ME = 0.0
    if Λ==Λ′ && F == F′
        ME = (-1)^(I+J′+F+S-Σ+J-Ω) * sqrt(I*(I+1)*(2I+1)*(2J+1)*(2J′+1)*S*(S+1)*(2S+1)) * 
            wigner6j_(J′,I,F,I,J,1) * 
            sum(
                wigner3j_(J,1,J′,-Ω,q,Ω′) * wigner3j_(S,1,S,-Σ,q,Σ′) 
                for q in -1:1
            )
    end

    return ME
end
export HyperfineIS



function TDM_E1(state::LinearCaseA, state′::LinearCaseA, p=nothing)
    """
    Transition dipole moment operator (basically Stark operator).
    Polarization has spherical component labeled by p.
    """

    Λ, S, Σ, J, Ω, I, F = unpack(state)
    Λ′, S′, Σ′, J′, Ω′, I′, F′ = unpack(state′)

    ME = 0.0
    if Σ==Σ′
        ME = (-1)^(J′+S+F+1) * wigner6j_(J′,F′,S,F,J,1) * sqrt((2F′+1)*(2F+1)) *
            (-1)^(J′-Ω′) * sqrt((2J′+1)*(2J+1)) * 
            sum(
                wigner3j_(J′,1,J,-Ω′,q,Ω)
                for q in -1:1
            )
    end

    return ME
end
export TDM_E1





