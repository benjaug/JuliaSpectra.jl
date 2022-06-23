using Parameters
using CompositeStructs
using InfiniteArrays
 
abstract type AbstractCaseB <: LinearBasisState end
export AbstractCaseB

export 
    LinearCaseB_Field,
    unpack,
    Rotation, SpinRotation, Hyperfine_IS, Hyperfine_Dipolar, ℓDoubling, Stark, Zeeman


Base.@kwdef struct LinearCaseB_Field <: AbstractCaseB
    Λ::Rational{Int64}
    N::Rational{Int64}
    S::Rational{Int64}
    J::Rational{Int64}
    I::Rational{Int64}
    F::Rational{Int64}
    M::Rational{Int64}

    constraints = (
        N = abs(Λ):∞,
        J = abs(N-S):abs(N+S),
        F = abs(J-I):abs(J+I),
        M = -F:F
    )
end


function unpack(state::LinearCaseB_Field)
    return (state.Λ, state.N, state.S, state.J, state.I, state.F, state.M)
end


function Rotation(state::LinearCaseB_Field, state′::LinearCaseB_Field)


    Λ, N, S, J, I, F, M = unpack(state)
    Λ′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = N*(N+1) - Λ^2

    return ME * δ(Λ,Λ′) * δ(N,N′) * δ(J,J′) * δ(F,F′) * δ(M,M′) * δ(S,S′) * δ(I,I′)
end


# Define the spherical tensor T^k_q(ϵ), here for linear and symmetric top molecules
const T = [
    0 0 0
    -1/√3 0 -1/√6
    0 0 0
    ]

function SpinRotation(state::LinearCaseB_Field, state′::LinearCaseB_Field)
    # See Hirota, eq. (2.3.35). There's a simpler form if Λ=0, but this will reduce to it when appropriate.
    Λ, N, S, J, I, F, M = unpack(state)
    Λ′, N′, S′, J′, I′, F′, M′ = unpack(state′)

    ME = (-1)^(J + S + N) * (-1)^(N - Λ) *
        sqrt( S * (S + 1) * (2S + 1) * (2N + 1) * (2N′ + 1) ) *
        wigner6j_(N′, S ,J, S, N, 1) *
        sum( sqrt(2k + 1) *
            (
                (-1)^k * 
                sqrt( N′ * (N′ + 1) * (2N′ + 1) ) * 
                wigner6j_(1, 1, k, N, N′, N′) +
                sqrt( N * (N + 1) * (2N + 1) ) *
                wigner6j_(1, 1, k, N′, N, N)
            ) *
                wigner3j_(N, k, N′, -Λ, q, Λ′) * T[q+2, k+1]
            for k in 0:2, q in -1:1
        )

    return  ME * δ(J, J′) * δ(F, F′) * δ(M, M′) * δ(S,S′) * δ(I,I′)
end


function Hyperfine_IS(state::LinearCaseB_Field, state′::LinearCaseB_Field)
    # Fermi-contact interaction
    # Hirota, pg. 39
    Λ, N, S, J, I, F, M = unpack(state)
    Λ′, N′, S′, J′, I′, F′, M′ = unpack(state′)
    return (-1)^(J′ + F + I + J + N + S + 1) *
        sqrt( (2J′ + 1) * (2J + 1) * S * (S + 1) * (2S + 1) * I * (I + 1) * (2I + 1) ) *
        wigner6j_(I, J, F, J′, I, 1) *
        wigner6j_(S, J, N, J′, S, 1) *
        δ(Λ, Λ′) * δ(N, N′) * δ(F, F′) * δ(M, M′)
end

function Hyperfine_Dipolar(state::LinearCaseB_Field, state′::LinearCaseB_Field)
    # Dipolar interaction term, from c(Iz ⋅ Sz)
    # Hirota, pg. 39
    Λ, N, S, J, I, F, M = unpack(state)
    Λ′, N′, S′, J′, I′, F′, M′ = unpack(state′)
    return - sqrt(30) * (-1)^(J′ + I + F + N) *
        wigner6j_(I, J, F, J′, I, 1) * 
        wigner9j_(N, N′, 2, S, S, 1, J, J′, 1) * 
        wigner3j_(N, 2, N′, -Λ, 0, Λ′) *
        sqrt( S * (S + 1) * (2S + 1) * I * (I + 1) * (2I + 1) * (2J + 1) * (2J′ + 1) * (2N + 1) * (2N′ + 1) ) *
        δ(F, F′) * δ(M, M′)
end


function ℓDoubling(state::LinearCaseB_Field, state′::LinearCaseB_Field)
    Λ, N, S, J, I, F, M = unpack(state)
    Λ′, N′, S′, J′, I′, F′, M′ = unpack(state′)
    return (-1)^(N - Λ) *
        (1 / (2 * sqrt(6))) *
        sqrt( (2N - 1) * (2N) * (2N + 1) * (2N + 2) * (2N + 3) ) *
        sum(
            wigner3j_(N, 2, N, -Λ, 2q, Λ′)
            for q in [-1,1]
        ) *
        δ(abs(Λ′ - Λ), 2) *
        δ(J, J′) * δ(F, F′) * δ(M, M′)
end


function Stark(state::LinearCaseB_Field, state′::LinearCaseB_Field)
    # Hirota, equation (2.5.35)
    Λ, N, S, J, I, F, M = unpack(state)
    Λ′, N′, S′, J′, I′, F′, M′ = unpack(state′)
    return -(-1)^(F - M) * 
        wigner3j_(F, 1, F′, -M, 0, M′) * 
        (-1)^(J + I + F′ + 1) * sqrt( (2F + 1) * (2F′ + 1) ) *
        wigner6j_(J, F, I, F′, J′, 1) *
        (-1)^(N + S + J′ + 1) * sqrt( (2J + 1) * (2J′ + 1) ) *
        wigner6j_(N, J, S, J′, N′, 1) *
        (-1)^(N - Λ) * sqrt( (2N + 1) * (2N′ + 1) ) *
        wigner3j_(N, 1, N′, -Λ, 0, Λ′)
end

function Zeeman(state::LinearCaseB_Field, state′::LinearCaseB_Field)
    # Hirota, equation (2.5.16) and (2.5.19)
    Λ, N, S, J, I, F, M = unpack(state)
    Λ′, N′, S′, J′, I′, F′, M′ = unpack(state′)
    return (-1)^(F - M) *
        wigner3j_(F, 1, F′, -M, 0, M) *
        (-1)^(J + I + F′ + 1) * sqrt( (2F + 1) * (2F′ + 1) ) * 
        wigner6j_(J′, F′, I, F, J, 1) *
        (-1)^(N + S + J + 1) * 
        sqrt( (2J + 1) * (2J′ + 1) * S * (S + 1) * (2S + 1) ) *
        wigner6j_(S, J′, N, J, S, 1) * 
        δ(Λ, Λ′) * δ(N, N′) * δ(M, M′)
end