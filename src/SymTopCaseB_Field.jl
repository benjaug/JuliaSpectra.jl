using Parameters
using CompositeStructs
using InfiniteArrays

""" 
Code for a C3v symmetric top molecule wtih S=1/2.
"""

Base.@kwdef struct SymTopCaseB_Field <: SymTopBasisState
    Λ::Rational{Int64}
    K::Rational{Int64}
    N::Rational{Int64}
    S::Rational{Int64}
    J::Rational{Int64}
    I::Rational{Int64}
    F::Rational{Int64}
    M::Rational{Int64}

    constraints = (
        N = 0:abs(K),
        J = abs(N-S):abs(N+S),
        F = abs(J-I):abs(J+I),
        M = -F:F
    )
end
export SymTopCaseB_Field

function unpack(state::SymTopCaseB_Field)
    return (state.Λ, state.K, state.N, state.S, state.J, state.I, state.F, state.M)
end
export unpack


function Rotation_K(state::SymTopCaseB_Field, state′::SymTopCaseB_Field)
    """
    K^2 part of the rotational Hamiltonian.
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = K^2

    return ME * δ(Λ,Λ′) * δ(N,N′) * δ(K,K′) * δ(F,F′) * δ(M,M′) * δ(J,J′) * δ(S,S′) * δ(I,I′)
end
export Rotation_K


function Rotation_N(state::SymTopCaseB_Field, state′::SymTopCaseB_Field)
    """
    N^2 part of the rotational Hamiltonian.
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = N*(N+1)

    return ME * δ(Λ,Λ′) * δ(N,N′) * δ(K,K′) * δ(F,F′) * δ(M,M′) * δ(J,J′) * δ(S,S′) * δ(I,I′)
end
export Rotation_N


function SpinRotation_1(state::SymTopCaseB_Field, state′::SymTopCaseB_Field)
    """
    Spin-rotation term that multiplies ϵbc
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = 1/2*(J*(J+1) - N*(N+1) - S*(S+1))

    return ME * δ(Λ,Λ′) * δ(N,N′) * δ(K,K′) * δ(F,F′) * δ(M,M′) * δ(J,J′) * δ(S,S′) * δ(I,I′)
end
export SpinRotation_1


function SpinRotation_2(state::SymTopCaseB_Field, state′::SymTopCaseB_Field)
    """
    Spin-rotation term that multiplies ϵaa - ϵbc
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = K*(2N+1)*sqrt(S*(S+1)*(2S+1)) * (-1)^(N+J+S) *
        wigner6j_(S,N,J,N,S,1) * wigner3j_(N,1,N,-K,0,K)

    return ME * δ(Λ,Λ′) * δ(N,N′) * δ(K,K′) * δ(F,F′) * δ(M,M′) * δ(J,J′) * δ(S,S′) * δ(I,I′)
end
export SpinRotation_2


function Hyperfine_Fermi(state::SymTopCaseB_Field, state′::SymTopCaseB_Field)
    """
    Fermi contact hyperfine term. Multiplies aF.
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = (-1)^(J+F+I) * wigner6j_(I,J,F′,J′,I,1) * sqrt(I*(I+1)*(2I+1)) *
        sqrt(S*(S+1)*(2S+1)) * (-1)^(J′+N′+1+S) * sqrt((2J′+1)*(2J+1)) *
        wigner6j_(S,J,N′,J′,S,1)

    return ME * δ(F,F′) * δ(M,M′) * δ(N,N′) * δ(K,K′) * δ(S,S′) * δ(I,I′) * δ(Λ,Λ′)
end
export Hyperfine_Fermi

function Hyperfine_Dipolar_TbbmTcc(state::SymTopCaseB_Field, state′::SymTopCaseB_Field)
    """
    Dipolar hyperfine term, multiplies (Tbb-Tcc)
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = (-1) * (-1)^(J+I+F′) * sqrt(30) * sqrt(S*(S+1)*(2S+1)*(2J′+1)*(2J+1)) *
        sqrt((2N′+1)*(2N+1)) * wigner6j_(I,J,F′,J′,I,1) * wigner9j_(N′,N,2,S,S,1,J′,J,1) *
        (-2*sqrt(I*(I+1)*(2I+1))) *
        sum(
            (-1)^(N′-K′) * 1/sqrt(24) * wigner3j_(N′,2,N,-K′,q,K)
            for q in [-2,2]
        )

    return ME * δ(F,F′) * δ(M,M′) * δ(S,S′) * δ(I,I′) * δ(Λ,Λ′)
end
export Hyperfine_Dipolar_TbbmTcc


function Hyperfine_Dipolar_Taa(state::SymTopCaseB_Field, state′::SymTopCaseB_Field)
    """
    Dipolar hyperfine term, multiplies Taa.
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = (-1) * (-1)^(J+I+F′) * sqrt(30) * sqrt(S*(S+1)*(2S+1)*(2J′+1)*(2J+1)) *
        sqrt((2N′+1)*(2N+1)) * wigner6j_(I,J′,F,J,I,1) * wigner9j_(N′,N,2,S,S,1,J′,J,1) *
        sqrt(I*(I+1)*(2I+1)) * (-1)^(N′-K′) * (1/2) * wigner3j_(N′,2,N,-K′,0,K)

    return ME * δ(F,F′) * δ(M,M′) * δ(S,S′) * δ(I,I′) * δ(Λ,Λ′)
end
export Hyperfine_Dipolar_Taa


function Stark(state::SymTopCaseB_Field, state′::SymTopCaseB_Field, p=0)
    """
    Stark matrix element. Assumes q=0 component (dipole along a-axis) and 
    takes an input p to specify direction of field relative to Z axis. Default
    to p=0.
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = (-1)^(p) * (-1)^(F′-M′) * wigner3j_(F′,1,F,-M′,p,M) *
        (-1)^(F+J′+1+I) * sqrt((2F+1)*(2F′+1)) * wigner6j_(J′,F′,I,F,J,1) * 
        (-1)^(J+N′+1+S) * sqrt((2J+1)*(2J′+1)) * wigner6j_(N′,J′,S,J,N,1) * 
        (-1)^(N′-K′) * sqrt((2N′+1)*(2N+1)) * wigner3j_(N′,1,N,-K′,0,K)


    return ME * δ(S,S′) * δ(I,I′) * δ(Λ,Λ′)
end
export Stark


function Zeeman_S(state::SymTopCaseB_Field, state′::SymTopCaseB_Field)
    """
    Zeeman matrix element for electron spin. Assume field along Z-axis.
    """

    Λ, K, N, S, J, I, F, M = unpack(state)
    Λ′, K′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = 2*(-1)^(N+S+J′+1) * sqrt((2J′+1)*(2J+1)*S*(S+1)*(2S+1)) * wigner6j_(S,J′,N,J,S,1) *
        (-1)^(F′-M) * wigner3j_(F′,1,F,-M,0,M) * (-1)^(J′+I+F+1) * sqrt((2F′+1)*(2F+1)) *
        wigner6j_(J′,F′,I,F,J,1)


    return ME * δ(K,K′) * δ(N,N′) * δ(M,M′)
end
export Zeeman_S